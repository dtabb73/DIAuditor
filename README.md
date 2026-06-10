# DIAuditor
Quality metrics generator for Data-Independent Acquisition proteomics

Quality metrics can play at least two key roles in Data-Independent Acquisition experiments.  First, the development and refinement of methods for DIA benefits substantially from having the ability to quantitative characterize the result of a new DIA method.  Second, when one has retrieved the raw data for a DIA experiment performed elsewhere, it can be a bit difficult to extract the method from those files.  DIAuditor is able to speed both of these goals substantially.  The software, created in C# and employing .NET libraries, can be run on the command line in Microsoft Windows (and hopefully I will soon have updated its GUI to work with the Visual Studio build system).  The software will process through all mzMLs in the current working directory.

DIAuditor reports metrics at two different levels: "DIAuditor-byRun.tsv" supplies one row of metrics for each input mzML, while "DIAuditor-byIsolationWindow.tsv" reports a row of metrics for each isolation window in each input mzML.  The metrics fall in these categories:

# Metrics by Run
+SourceFile: the root name of the LC-MS/MS experiment (such as RAW or directory name, without extension)
+Instrument: the model of instrument reported in the mzML file
+SerialNumber: the serial number of the instrument, if reported in the mzML file
+StartTimeStamp: the time and date when this LC-MS/MS began acquisition
+RTDuration: the number of minutes into the run when scan acquisition stopped
+mzMLMS1Count: the number of MS scans recorded for the run
+mzMLMSnCount: the number of MSn scans recorded for the run
+MS1Resolution: the mass resolving power recorded for the last MS scan
+MS1TIC25ileRT, MS1TIC50ileRT, MS1TIC75ileRT: the retention times at which the first, second, and third quartiles of the TIC (total ion current) sum across MS scans are accumulated
+MS1TotalTIC: the sum of TIC values for all MS scans in the run
+MS1CycleTime: the median time in seconds between successive MS scans in the run
+MS1PkCountMin, MS1PkCount25ile, MS1PkCount50ile, MS1PkCount75ile, MS1PkCountMax: the quartiles from the distribution of centroid peak counts from all MS scans in the run
+IsolationWindowCount: the number of isolation windows in the method, requiring minimum and maximum m/z boundaries as well as ion mobility to be identical
+CyclesMin and CyclesMax: the range in the count of MSn scans for all isolation windows
+MZRangeMin and MZRangeMax: the minimum m/z for isolation window definitions and the maximum m/z for isolation window definitions
+IsolationWindowWidthMin and IsolationWindowWidthMax: the minimum and maximum width in m/z for any isolation window
+AverageMedianCycleTime: the average for all isolation window median cycle times
+TICMedianRTMin and TICMedianRTMax: the range of retention times at which half of all isolation window TIC has been accumulated
+TotalTiCMin and TotalTICMax: the range of TICs accumulated across all MSn scans of each isolation window
+PkCountMedianMin and PkCountMedianMax: the range of median peak counts for all MSn scans of each isolation window
# Metrics by Isolation Window
+LoMZ, HiMZ: the m/z bounds of this isolation window
+WidthMZ: the size of this isolation window (HiMZ - LoMZ)
+IonMobility: the FAIMs energy for this isolation window
+MassResolvingPower: the resolution at which the MSns of this isolation window were acquired
+MSMSCount: the number of MSns acquired for this isolation window
+RTMin and RTMax: the range of retention times during which this isolation window was measured
+CycleTimeMedian: the median number of seconds between successive MSn scans in this isolation window
+TIC25ileRT, TIC50ileRT, TIC75ileRT: the retention times at which the first, second, and third quartiles of the TIC sum across MSn scans for this isolation window are accumulated
+TotalTIC: the sum of TIC for all MSn scans in this isolation window
+PkCountMin, PkCount25ile, PkCount50ile, PkCount75ile, PkCountMax: the quartiles from the distribution of centroid peak counts from all MSn scans in this isolation window
