from pyproj import Geod
import numpy as np
import pandas as pd
from matplotlib.path import Path
import streamlit as st
import io

def obstacle_check(icao, rwy, obs_lat, obs_lon, height, df):
    
    def destination_point(lat0, lon0, bearing_deg, distance_m):
        """
        Calculates the destination point from an initial geodetic position,
        given a bearing (degrees) and distance (meters).
    
        Returns:
        - lat1, lon1: destination coordinates (degrees)
        """
        geod = Geod(ellps="WGS84")
        lon1, lat1, _ = geod.fwd(lon0, lat0, bearing_deg, distance_m)
        return lat1, lon1
    
    def get_airport_info(df, icao, rwy):
      thr_lon = df.loc[(df['Icao'] == icao) & (df['Name4'] == rwy) & (df['Ident'].isnull())]['ThresholdLongitude'].values[0]
      thr_lat = df.loc[(df['Icao'] == icao) & (df['Name4'] == rwy) & (df['Ident'].isnull())]['ThresholdLatitude'].values[0]
      magnetic_variation = df.loc[(df['Icao'] == icao)]['MagneticVariation'].values[0]
      direction = magnetic_variation[-1]
      magnetic_variation = int(magnetic_variation[:-1])/10000
      if direction in ['S', 'W']:
        magnetic_variation *= -1
      asda = df.loc[(df['Icao'] == icao) & (df['Name4'] == rwy) & (df['Ident'].isnull())]['Asda'].values[0]
      toda = df.loc[(df['Icao'] == icao) & (df['Name4'] == rwy) & (df['Ident'].isnull())]['Toda'].values[0]
      tora = df.loc[(df['Icao'] == icao) & (df['Name4'] == rwy) & (df['Ident'].isnull())]['Tora'].values[0]
      lda = df.loc[(df['Icao'] == icao) & (df['Name4'] == rwy) & (df['Ident'].isnull())]['Lda'].values[0]
      magnetic_heading = df.loc[(df['Icao'] == icao) & (df['Name4'] == rwy)]['MagneticHeading'].values[0]
      if magnetic_heading > int(rwy[1:2])*10 + 20 or magnetic_heading < int(rwy[1:2])*10 - 20:
        magnetic_heading = (magnetic_heading + 180) % 360
      thr_elevation = df.loc[(df['Icao'] == icao) & (df['Name4'] == rwy)]['ThresholdElevation'].values[0]
      thr_elevation = thr_elevation*3.28084
      true_heading = magnetic_heading + magnetic_variation

      def convert_to_decimal(coord: str) -> float:
        coord = coord.strip()
        direction = coord[-1]
        coord = coord[:-1]

        if len(coord) == 8:
            deg = int(coord[0:2])
            minute = int(coord[2:4])
            sec = int(coord[4:])
        elif len(coord) == 9:
            deg = int(coord[0:3])
            minute = int(coord[3:5])
            sec = int(coord[5:])
        elif len(coord) == 6:
            deg = 0
            minute = int(coord[0:2])
            sec = int(coord[2:])
        elif len(coord) == 7:
            deg = int(coord[0])
            minute = int(coord[1:3])
            sec = int(coord[3:])
        else:
            raise ValueError("Error Coordinate")

        decimal_deg = deg + (minute / 60) + (sec / 360000)

        if direction in ["S", "W"]:
            decimal_deg = -decimal_deg

        return decimal_deg

      thr_lon = convert_to_decimal(thr_lon)
      thr_lat = convert_to_decimal(thr_lat)

      end_of_tora_lat, end_of_tora_lon = destination_point(thr_lat, thr_lon, true_heading, lda)

      return thr_lon, thr_lat, asda, toda, tora, lda, true_heading, end_of_tora_lat, end_of_tora_lon, thr_elevation

    def funnel_check(thr_lat, thr_lon, toda, obs_lat, obs_lon, true_heading):

        x0, y0 = destination_point(thr_lat, thr_lon, true_heading, 500)
        x1, y1 = destination_point(x0, y0, (true_heading - 90) % 360, 90)
        x2, y2 = destination_point(x0, y0, (true_heading + 90) % 360, 90)
        x3, y3 = destination_point(x1, y1, true_heading, toda-500)
        x4, y4 = destination_point(x2, y2, true_heading, toda-500)
        x5, y5 = destination_point(x3, y3, (true_heading - 7.1) % 360, 1699)
        x6, y6 = destination_point(x4, y4, (true_heading + 7.1) % 360, 1699)
        x7, y7 = destination_point(x5, y5, true_heading, 100000)
        x8, y8 = destination_point(x6, y6, true_heading, 100000)

        shapeX = [x0, x1, x2, x3, x4, x5, x6, x7, x8]
        shapeY = [y0, y1, y2, y3, y4, y5, y6, y7, y8]

        def inpolygon(x, y, poly_x, poly_y):
          points = np.vstack((x, y)).T
          polygon = np.vstack((poly_x, poly_y)).T
          path = Path(polygon)
          return path.contains_points(points)

        inside = inpolygon([obs_lat], [obs_lon], shapeX, shapeY)[0]
        return inside

    def obstacle_distances(obs_lat, obs_lon, thr_lat, thr_lon, end_of_tora_lat, end_of_tora_lon, true_heading):

      def distance_and_bearing(lat1, lon1, lat2, lon2):
        """
        Calculates the distance (in meters) and initial bearing (degrees from North)
        from point (lat1, lon1) to point (lat2, lon2).
        """
        geod = Geod(ellps="WGS84")
        az12, az21, dist = geod.inv(lon1, lat1, lon2, lat2)
        return dist, az12

      distance_from_start_of_tora, bearing_from_start_of_tora = distance_and_bearing(thr_lat, thr_lon, obs_lat, obs_lon)
      distance_from_end_of_tora, bearing_from_end_of_tora = distance_and_bearing(end_of_tora_lat, end_of_tora_lon, obs_lat, obs_lon)
      bearing_from_start_of_tora = (bearing_from_start_of_tora + 360) % 360
      bearing_from_end_of_tora = (bearing_from_end_of_tora + 360) % 360

      longitudinal_distance_start_of_tora = distance_from_start_of_tora * np.cos(np.deg2rad(abs((-bearing_from_start_of_tora + true_heading) % 360)))
      longitudinal_distance_end_of_tora = distance_from_end_of_tora * np.cos(np.deg2rad(abs((-bearing_from_end_of_tora + true_heading) % 360)))

      return longitudinal_distance_start_of_tora, longitudinal_distance_end_of_tora
    thr_lon, thr_lat, asda, toda, tora, lda, true_heading, end_of_tora_lat, end_of_tora_lon, thr_elevation = get_airport_info(df, icao, rwy)
    height_from_threshold = height - thr_elevation
    inside = funnel_check(thr_lat, thr_lon, toda, obs_lat, obs_lon, true_heading)
    longitudinal_distance_start_of_tora, longitudinal_distance_end_of_tora = obstacle_distances(obs_lat, obs_lon, thr_lat, thr_lon, end_of_tora_lat, end_of_tora_lon, true_heading)
    return inside, longitudinal_distance_start_of_tora, longitudinal_distance_end_of_tora, height_from_threshold


st.title('Aerodrome Obstacle Analysis Tool ✈️')

DATABASE_PATH = 'obstacle_distance/data/TEST.xlsx'
# Input fields
try:
    df_database = pd.read_excel(DATABASE_PATH, engine='openpyxl')
    st.success('Database loaded successfully!')

    # Input fields
    icao = st.text_input('ICAO Code', 'EDDM')
    rwy = st.text_input('Runway Name', '08L')
    obs_lat = st.number_input('Obstacle Latitude', format="%.6f")
    obs_lon = st.number_input('Obstacle Longitude', format="%.6f")
    height = st.number_input('Obstacle Height (ft)')

    if st.button('Run Analysis'):
        is_inside, dist_start, dist_end, height_from_threshold = obstacle_check(
            icao, rwy, obs_lat, obs_lon, height, df_database
        )
        
        # Display the results
        if is_inside:
            st.success('✅ The obstacle is INSIDE the takeoff funnel.')
        else:
            st.warning('⚠️ The obstacle is OUTSIDE the takeoff funnel.')
            
        st.info(f'**Longitudinal distance from runway start:** {dist_start:.2f} meters')
        st.info(f'**Longitudinal distance from TORA end:** {dist_end:.2f} meters')

except FileNotFoundError:
    st.error(f"Database file not found. Please ensure '{DATABASE_PATH}' is in your GitHub repository.")
except Exception as e:
    st.error(f"An error occurred: {e}")
