# DIAuditor
Quality metrics generator for Data-Independent Acquisition proteomics

Quality metrics can play at least two key roles in Data-Independent Acquisition experiments.  First, the development and refinement of methods for DIA benefits substantially from having the ability to quantitative characterize the result of a new DIA method.  Second, when one has retrieved the raw data for a DIA experiment performed elsewhere, it can be a bit difficult to extract the method from those files.  DIAuditor is able to speed both of these goals substantially.  The software, created in C# and employing .NET libraries, can be run on the command line in Microsoft Windows (and hopefully I will soon have updated its GUI to work with the Visual Studio build system).  The software will process through all mzMLs in the current working directory.

DIAuditor reports metrics at two different levels: "DIAuditor-byRun.tsv" supplies one row of metrics for each input mzML, while "DIAuditor-byIsolationWindow.tsv" reports a row of metrics for each isolation window in each input mzML.  The metrics fall in these categories:
+ Counts and Reports: Features like the number of MS and MS/MS scans or the "mass resolving power" of scans will be reported along with identifying information like the instrument model and serial number.
+ Window definitions: The lower bound m/z, target m/z, and upper bound m/z as well as FAIMS CV for each isolation window is reported, and many metrics specific to each window are reported.
+ Retention Time and quartiles of TIC: We frequently see that different isolation windows will vary in the retention time at which the TIC maximizes; high m/z windows are likely to see their TIC values crest later in retention time.
+ Cycle Times: the rate at which MS scans are acquired may differ from the rate at which MS/MS scans are acquired in a given isolation window.  The median time in seconds is reported at both levels.
+ Peak Count Quartiles: Employing features like FAIMS can have a substantial impact on how many peaks are centroided from each MS or MS/MS scan.  These are reported for each isolation window separately.
