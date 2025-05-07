# -*- coding: utf-8 -*-
"""
Pia Bhatia 
Drndic Lab 
11-5-2024 

Created with help of OpenAI ChatGPT. This code opens up a TIF file, allows you to
click in the center of a pore, then it will draw 4 line profiles with a specific 
integration width (in units of pixels). It then produces a plot of the normalized 
intensities along this line profile. Both figures are saved as svgs and pngs in the 
specified directories. The code also draws vertical lines where the image intensity 
in each profile crosses a specific threshold value. The code computes the distance between 
the two vertical lines, and averages them to report a pore diameter with error (stdev). 
"""

# -*- coding: utf-8 -*-

# Required imports
import rasterio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection  # Import for line width control
import datetime

# Close all open matplotlib figures
plt.close('all')

# Define the path and parameters
path = r'C:/Users/19085/OneDrive/Desktop/GURU 2.0/Figures/Device J #/ADF_200kV_20MX_GURU2.0_CTRL_ChipJ_0030.tif'
savedirectory = path 
pixelsize = 10/1024  # in units of nm per pixel
intensitythresh = 0.2  # Threshold between 0 and 1 (normalized); 20% of intensity 
integrationwidth = 5 
linewidth = 5  # in units of pixels (same as integration width)
linelength = 4 #length of line profile in units of nm 
SBSize = 3 #in units of nm 
letter = 'J'

#To draw the lines on the image with the correct thickness, we have to perform the following conversion. 
#In the 'oneclick(event)'function we use LineCollection() to draw the lines on the image. 
#LineCollection() has a parameter called 'linewidth' and the default unit for this parameter is in POINTS. 
#POINTS are NOT the same as pixels. Given that the figure DPI = 100, we know that there are 100 dots (pixels) per inch.
#1 point = 1/72 inch  
# points = (pixels / DPI) * 72. 
#This gives us the following conversion: 
integrationwidth_units_points = (linewidth / 100) * 72 #for DPI = 100 
#You can verify this conversion but opening the .svg file in adobe illustrator and zooming in on a line profile. 
#You'll see that it is the desired number of pixels wide. 

color1 = '#481567FF'
color2 = '#EFCC00'
color3 = '#287D8EFF'
color4 = '#8ACE00'  # brat green hehe
ImageWithLinesNameSVG = f"Device{letter}_WithLineProfilesDrawn_{SBSize}nmSB.svg"
ImageWithLinesNamePNG = f"Device{letter}_WithLineProfilesDrawn_{SBSize}nmSB.png"
ImageWithLinesNameTIF = f"Device{letter}_WithLineProfilesDrawn_{SBSize}nmSB.tif"
LineProfilePlotNameSVG = f"Device{letter}_LineProfile_4nmLongLines_{intensitythresh}_IntensityThreshold_{linewidth}px_IntegrationWidth.svg"
LineProfilePlotNamePNG = f"Device{letter}_LineProfile_4nmLongLines_{intensitythresh}_IntensityThreshold_{linewidth}px_IntegrationWidth.png"
LineProfilePlotNameTIF = f"Device{letter}_LineProfile_4nmLongLines_{intensitythresh}_IntensityThreshold_{linewidth}px_IntegrationWidth.tif"

# Set plot font properties (Arial, size 12, and bold axis titles)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16

# Function to create a custom colormap from specified colors
def generate_custom_colormap(colors):
    return mcolors.LinearSegmentedColormap.from_list("custom_gradient", colors, N=len(colors))

# Function to draw line profiles and calculate average intensity with integration width
def draw_line_profile(image, start_point, end_point, integration_width=integrationwidth):
    """
    Draws the profile of pixel intensities along a line, following the true diagonal direction.
    The integration width defines how wide the sampling window is along the line.
    """
    y_start, x_start = start_point
    y_end, x_end = end_point
    length = np.linalg.norm(np.array(end_point) - np.array(start_point))
    
    profile = []
    num_samples = int(length)  # Number of points to sample along the line

    # Calculate the direction vector (dx, dy)
    dx = x_end - x_start
    dy = y_end - y_start

    # Normalize the direction vector (dx, dy)
    norm = np.sqrt(dx**2 + dy**2)
    dx /= norm
    dy /= norm

    for i in range(num_samples):
        # Move along the line by i steps
        x = int(x_start + i * dx)
        y = int(y_start + i * dy)

        # Check if the point is within the image boundaries
        if 0 <= x < image.shape[1] and 0 <= y < image.shape[0]:
            # Define the window around the point (integration_width applies in both directions)
            window_start_x = max(0, x - integration_width // 2)
            window_end_x = min(image.shape[1], x + integration_width // 2)
            window_start_y = max(0, y - integration_width // 2)
            window_end_y = min(image.shape[0], y + integration_width // 2)

            # Extract the intensity values in the window
            intensity_window = image[window_start_y:window_end_y, window_start_x:window_end_x]

            # Calculate the average intensity in the window
            intensity = np.mean(intensity_window)
            profile.append(intensity)

    return np.array(profile)

# Load the TIFF file
def load_tif(file_path):
    with rasterio.open(file_path) as src:
        image = src.read(1)
        return image

# Convert length in nanometers to pixels
def nanometers_to_pixels(length_nm, pixel_size_nm):
    return int(length_nm / pixel_size_nm)

# Save the modified image with drawn lines
def save_image_with_lines(original_image, lines, file_path):
    output_image = original_image.copy()
    for line in lines:
        start, end, color = line
        y_start, x_start = start
        y_end, x_end = end
        rr, cc = np.linspace(y_start, y_end, num=100).astype(int), np.linspace(x_start, x_end, num=100).astype(int)
        output_image[rr, cc] = np.clip(output_image[rr, cc] + color[0], 0, 255)

    with rasterio.open(file_path, 'w', driver='GTiff', height=output_image.shape[0], width=output_image.shape[1], count=1, dtype='uint8') as dst:
        dst.write(output_image, 1)


#function to find the crossing points based on NORMALIZED INTENSITY of LP 
def find_crossing_points(profile, threshold, midpoint_index): 
    # Initialize left_crossing and right_crossing to None in case no crossings are found
    left_crossing = None
    right_crossing = None
    normprofile = (profile - np.min(profile))/(np.max(profile)-np.min(profile))
    
    # Search for the leftmost crossing point before the midpoint
    for i in range(midpoint_index - 1, -1, -1):  # Search leftwards (before midpoint)
        if normprofile[i] < threshold and normprofile[i-1] >= threshold:
            # Interpolate the exact crossing point
            left_crossing = i + (threshold - normprofile[i]) / (normprofile[i-1] - normprofile[i])
            break  # Once the first crossing is found, exit loop

    # Search for the rightmost crossing point after the midpoint
    for i in range(midpoint_index, len(normprofile) - 1):  # Search rightwards (after midpoint)
        if normprofile[i] < threshold and normprofile[i+1] >= threshold:
            # Interpolate the exact crossing point
            right_crossing = i + (threshold - normprofile[i]) / (normprofile[i+1] - normprofile[i])
    
    # Return the leftmost and rightmost crossing points
    return left_crossing, right_crossing
    
# Main code
if __name__ == "__main__":
    file_path = path
    pixel_size_nm = pixelsize
    line_length_nm = linelength 
    integration_width = linewidth
    intensity_threshold = intensitythresh
    line_length_pixels = nanometers_to_pixels(line_length_nm, pixel_size_nm)
    angles = [0, 45, 90, 135]

    # Create the custom colors array and colormap
    custom_colors = [color1, color2, color3, color4]  # Your custom colors
    colormap = generate_custom_colormap(custom_colors)

    original_image = load_tif(file_path)
    
    fig, ax = plt.subplots(figsize=(10,10), dpi=100) #figsize is in units of inches 

    # Remove axes' ticks, labels, and borders
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Adjust the layout to ensure the image fills the entire window
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    ax.imshow(original_image, cmap='gray')
    # Set aspect ratio to 'equal' to ensure pixel-based linewidths are correct
    ax.set_aspect('equal', adjustable='box')
    plt.axis('off')

    # Optionally set y-limits (0 to image height)
    #for some reason without these lines axes went from -0.5 to 1023.5;
    #this corrects that! 
    plt.ylim(original_image.shape[0], 0) 
    plt.xlim(0, original_image.shape[1])

    state = {
        'start_point': None,
        'end_point': None,
        'profiles': [],
        'lines': []
    }

    distances = []  # List to store distances between crossing points

# Global variables to store information for text file
clicked_coordinates = None  # Stores x, y from the click event
average_distance = None     # Stores the calculated average distance
std_distance = None         # Stores the standard deviation of the distances

def onclick(event):
    global clicked_coordinates, average_distance, std_distance  # Declare global variables

    if event.xdata is not None and event.ydata is not None:
        x = int(event.xdata)
        y = int(event.ydata)
        
        clicked_coordinates = (x, y)  # Store clicked coordinates
        
        # Print the x and y coordinates
        print(f"Mouse clicked at (x, y): ({x}, {y})")
        
        if state['start_point'] is None:
            state['start_point'] = (y, x)
            ax.plot(x, y, 'ro')  # Plot a red dot at the point selected
            
        else:
            midpoint = state['start_point']  # Midpoint of line is set to the point user clicks (start point)
            state['profiles'] = []
            state['lines'] = []
            distances.clear()  # Clear previous distances
             
            # For loop runs over each angle in angles 
            for i, angle in enumerate(angles):
                radians = np.radians(angle)  # Converts angle from degrees to radians
                new_end_x = midpoint[1] + line_length_pixels * np.cos(radians)
                new_end_y = midpoint[0] + line_length_pixels * np.sin(radians)
                new_start_x = midpoint[1] - line_length_pixels * np.cos(radians)
                new_start_y = midpoint[0] - line_length_pixels * np.sin(radians)

                line_color = colormap(i)  # Color from colormap already in [0, 1] range

                # Define the line segments using a 2D array of points for the linewidth
                line_segments = [[(new_start_x, new_start_y), (new_end_x, new_end_y)]]

                # Create a LineCollection with pixel-wide lines
                line_collection = LineCollection(line_segments, linewidths=integrationwidth_units_points, colors=[line_color])
                ax.add_collection(line_collection)
                
                profile = draw_line_profile(original_image, (new_start_y, new_start_x), (new_end_y, new_end_x), integration_width)
                state['profiles'].append((angle, profile, line_color))
                normprofile = (profile - np.min(profile))/(np.max(profile)-np.min(profile))

                # Find the left and right crossing points based on the threshold
                midpoint_index = len(normprofile) // 2  # Find the midpoint of the profile
                left_crossing, right_crossing = find_crossing_points(normprofile, intensity_threshold, midpoint_index)

                if left_crossing is not None and right_crossing is not None:
                    # Convert the index to pixel position (nm)
                    left_position = left_crossing * pixel_size_nm
                    right_position = right_crossing * pixel_size_nm

                    # Calculate distance between crossings in nm
                    distance = abs(right_position - left_position)
                    distances.append(distance)  # Store distance

            # Calculate average and standard deviation of distances
            if distances:
                average_distance = np.mean(distances)
                std_distance = np.std(distances)
                print(f'Distances between crossings: {distances}')
                print(f'Average distance: {average_distance:.2f} nm')
                print(f'Standard deviation of distances: {std_distance:.2f} nm')

            state['start_point'] = None
            plt.draw()
            plt.savefig(savedirectory + ImageWithLinesNameSVG)
            plt.savefig(savedirectory + ImageWithLinesNamePNG)
            plt.savefig(savedirectory + ImageWithLinesNameTIF)

            # Plot intensity profiles with vertical lines at crossing points
            plt.figure(figsize=(10, 5))
            for angle, profile, color in state['profiles']:
                normprofile2 = (profile-np.min(profile))/(np.max(profile)-np.min(profile))
                x_nm = np.arange(len(normprofile2)) * pixel_size_nm
                plt.plot(x_nm, normprofile2, label=f'LP{state["profiles"].index((angle, profile, color)) + 1}', color=color, linewidth=2)

                # Draw vertical lines at the crossing points
                midpoint_index = len(normprofile2) // 2  # Midpoint index of profile
                left_crossing, right_crossing = find_crossing_points(normprofile2, intensity_threshold, midpoint_index)

                if left_crossing is not None:
                    plt.axvline(x=left_crossing * pixel_size_nm, color=color, linestyle='--', linewidth=1)
                if right_crossing is not None:
                    plt.axvline(x=right_crossing * pixel_size_nm, color=color, linestyle='--', linewidth=1)

            plt.xlabel('Distance (nm)')
            plt.ylabel('Normalized Pixel Intensity')
            plt.legend()
            plt.grid()
            plt.show()
            plt.savefig(savedirectory + LineProfilePlotNameSVG)
            plt.savefig(savedirectory + LineProfilePlotNamePNG)
            plt.savefig(savedirectory + LineProfilePlotNameTIF)

            # Write the text file with details
            write_text_file()

def write_text_file():
    # Get the current date and time for file naming
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Prepare the content for the text file
    text_content = f"""
    GURU 2.0 DEVICE {letter} PARAMETERS 
    Date processed: {current_date}
    File processed: {path}
    Pixel size in image: {pixelsize} nm/px
    Intensity threshold for defining pore size: {intensitythresh}
    Line width/integration width used for plotting profiles: {linewidth} px
    Line length for profiles: {linelength} nm
    Coordinates of the clicked point: x = {clicked_coordinates[0]}, y = {clicked_coordinates[1]}
    Average pore diameter: {average_distance:.2f} nm
    Standard deviation of pore diameter: {std_distance:.2f} nm
    """

    # Define the text file name with date and time
    text_file_name = f"{savedirectory}Device{letter}_parameters_{current_date.replace(':', '_').replace(' ', '_')}.txt"
    
    # Write the content to the text file
    with open(text_file_name, 'w') as f:
        f.write(text_content)
        
# Connect the callback function to the click event
cid = fig.canvas.mpl_connect('button_press_event', onclick)

# Show the interactive plot
plt.show()
