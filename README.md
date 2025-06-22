## Data Extraction - Main_Extract_Data.m
This script is the primary pipeline for converting Medtronic JSON files into a structured MATLAB-based **DBS (Deep Brain Stimulation) data structure**.

### 📌 What It Does
- Loads or initializes a `DBS_data.mat` structure to store patient data.
- Recursively scans a selected folder (or processes a single file) to find all `.json` files.
- Extracts most types of data stored in the JSON files:
    - **TrendLogs** (home sensing)
    - **Events**
    - **Groups**
    -  **info**
    - **Streaming** 
    - **Montages**
- Organizes extracted information into per-patient substructures.
---
#### 📂 Folder and Data Structure
- Creates or updates a `DBS_data.mat` file in the chosen folder.
- Organizes data in a nested struct - DBS_data.Patient_XX.(Data Type).
* The extracted data from all JSON files belonging to a patient is organized as a 44 rows × _N_ columns table.

|       | 01 Jan 2025 | 02 Jan 2025 | .... | 01 July 2025 |
| ----- | ----------- | ----------- | ---- | ------------ |
| 00:00 |             |             |      |              |
| 00:10 |             |             |      |              |
| 00:20 |             |             |      |              |
| ..... |             |             |      |              |
| 23:50 |             |             |      |              |
* The 10 minute increments are the result of the device's LFP home sensing temporal resolution. A round down from the exact measurement time is made to create the table.
* This way, the data is continuous and comparable across types, patients, hemispheres.
* Similar tables are also created for:
  - Stimulation Amplitude
  - Stimulation Pulse Width
  - Stimulation Contact choice
  - Sensing Frequency

---

## Data Plotting - Main_Plot_Figures.m
This script provides a modular interface for **visualizing neural and stimulation data** extracted into `DBS_data.mat`. In it we use a set of plotting functions to explore the data.

---

###  Available Plotting Functions

#### 1. `Diurnal_LFP_Plot`
**Visualizes circadian cycles** in LFP activity for selected frequency bands using 24-hour polar plots.

**Inputs:**  
`DBS_data`, patient(s), hemisphere, frequency band, day range  
**Options:** Plot mean across patients, custom colors

>**Example:**
`Diurnal_LFP_Plot(DBS_data, {'Patient_1', 'Patient_2'}, 'Left', [13 30], {[1 10], [20 30]}, 'PlotMean', true);`

![image](https://github.com/user-attachments/assets/a5cb0407-cbed-4e72-985f-43ed199e9417)

#### 2. `Heatmap_Plot`
**Creates a time × day heatmap** of LFP power for a specific patient and hemisphere.

**Inputs:**  
`DBS_data`, patient name, hemisphere  
**Options:** Day range, custom color map

>**Example:**
`Heatmap_Plot(DBS_data, 'Patient_1', 'Right', [1 20]);`

#### 3. `Histogram_Plot`
**Plots distributions** of LFP amplitudes:
- Left vs Right
- Day vs Night

**Inputs:**  
`DBS_data`, patient name, frequency band  
**Options:** Day range (or different ranges per hemisphere)

>**Example:**
`Histogram_Plot(DBS_data, 'Patient_1', [13 30], [10 40]);'

![image](https://github.com/user-attachments/assets/8ab84a69-1a1f-4280-9633-b63ca182ce7c)


#### 4. `LFP_All_Plots`
**Plot LFP over time**. Automatically divides the data based on unique combinations of stimulation parameters.
- LFP amplitude
- Stimulation amplitude
- Contact and frequency configurations

**Inputs:**  
`DBS_data`, patient name(s), hemisphere

> Example:
`LFP_All_Plots(DBS_data, 'Patient_1', 'Right');`

![image](https://github.com/user-attachments/assets/2387770c-4ca0-4d61-a40a-4b629ac22507)

#### 5. `LFP_Daily_Plot`
**Plots daily mean time-series of LFP or TEED with optional overlays:
- Stimulation changes
- Sensing configs
- Clinical events

**Inputs:**  
`DBS_data`, patient name(s), hemisphere, frequency band, day range  
**Options:** Show TEED instead of LFP, show events by type or frequency band

> Example:
`LFP_Daily_Plot(DBS_data, 'Patient_1', 'Left', [13 30], [1 60], 'Show_Events', true, 'Plot_TEED', false);`

![image](https://github.com/user-attachments/assets/98d777fe-3ab7-4d8f-9aeb-082b8f1bdbef)
---

#### 6. `Periodogram_Plot`
**Generates PSD plots** across multiple days to reveal underlying oscillatory patterns.

**Inputs:**  
`DBS_data`, patient name, hemisphere, day range  
**Usage:** Understand longer-term cycles (e.g., circadian, ultradian)

> Example:
`Periodogram_Plot(DBS_data, 'Patient_1', 'Left', [1 60]);`

![image](https://github.com/user-attachments/assets/23303be8-bbd1-4ad2-aa85-2537da8ca378)


---
## 🚀 How to Run
You can run the script interactively in MATLAB. It will prompt you for:
- A target folder (if `use_current_path == 0`, otherwise the file's folder will be used)
- A folder or file selection depending on `Process_all`:
  - If `Process_all=false` : You will be prompted to select only one JSON file you would like to process.
  - If `Process_all=false` : You will be prompted to select a folder. All JSON files in that folder, or in it's sub folders will be processed.
    
To begin the process open `Main_Extract_Data.m`, select the options that fit your needs and run the file. 
This will process the files, print progress messages, and update or create the DBS structure file.

For plotting:
1. Ensure `DBS_data.mat` is present in your working directory.
2. Open `Main_Plot_Figures.m`.
3. Modify the script to call the plot functions you want to run.
4. Run the script or individual sections from MATLAB.
 - All plotting functions can be ran manually, they are in `Plot Functions/`
