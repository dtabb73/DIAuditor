using MzqcCsLib;
using System;
using System.ComponentModel;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Threading;
using System.Xml;

namespace DIAuditor
{
    static class Program
    {
        static void Main(string[] args)
        {
            const string Version = "20260415 beta";
            // Use periods to separate decimals
            CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
            Console.WriteLine("DIAuditor: Quality metrics for Data-Independent Acquisition experiments");
            Console.WriteLine("David L. Tabb, ERIBA of University Medical Center of Groningen");
            Console.WriteLine("Version " + Version);
            //Console.WriteLine("--DIANN: read report.pr_matrix.tsv peptide matrix");

            var ReadDIANN = false;
            foreach (var item in args)
            {
                switch (item)
                {
                    case "--DIANN":
                        ReadDIANN = true;
                        break;
                    default:
                        Console.Error.WriteLine("\tError: I don't understand this argument: {0}.", item);
                        break;
                }
            }
            var CWD = Directory.GetCurrentDirectory();
            const string mzMLPattern = "*.mzML";
            var mzMLs = Directory.GetFiles(CWD, mzMLPattern);
            var Raws = new LCMSMSExperiment();
            var RawsRunner = Raws;
            string Basename;
            Stopwatch Timer = new Stopwatch();
            TimeSpan Duration;
            Timer.Start();
            Console.WriteLine("\nImporting from mzML files...");
            foreach (var current in mzMLs)
            {
                Basename = Path.GetFileNameWithoutExtension(current);
                Console.WriteLine("\tReading mzML {0}", Basename);
                RawsRunner.Next = new LCMSMSExperiment();
                var FileSpec = Path.Combine(CWD, current);
                var XMLfile = XmlReader.Create(FileSpec);
                RawsRunner = RawsRunner.Next;
                RawsRunner.SourceFile = Basename;
                RawsRunner.FileURI = new Uri(current);
                RawsRunner.ReadFromMZML(XMLfile);
                RawsRunner.ParseScanNumbers();
            }
            Timer.Stop();
            Duration = Timer.Elapsed;
            Console.WriteLine("\tTime for mzML reading: {0}", Duration.ToString());
            if (ReadDIANN == true)
            {
                Timer.Reset();
                Timer.Start();
                Console.WriteLine("\nReading DIANN result is not implemented yet.");
                Timer.Stop();
                Duration = Timer.Elapsed;
                Console.WriteLine("\tTime for DIANN reading: {0}", Duration.ToString());
            }
            Timer.Reset();
            Timer.Start();
            Console.WriteLine("\nComputing metrics for each isolation window of each experiment...");
            RawsRunner = Raws.Next;
            while (RawsRunner != null)
            {
                Console.WriteLine("\t\t" + RawsRunner.SourceFile);
                RawsRunner.ComputeMS1Metrics();
                RawsRunner.ComputeMetricsForSwaths();
                RawsRunner = RawsRunner.Next;
            }
            Timer.Stop();
            Duration = Timer.Elapsed;
            Console.WriteLine("\tTime for isolation window evaluation: {0}", Duration.ToString());

            Console.WriteLine("\nWriting DIAuditor TSV reports...");
            Raws.WriteTextQCReport();
            Console.WriteLine("Writing DIAuditor.mzQC...");
            Raws.WriteMZQCReport(Version);
            Timer.Stop();
            Duration = Timer.Elapsed;
            Console.WriteLine("\tTime for reporting: {0}", Duration.ToString());
        }
    }

    class MS1Scan
    {
        public float TIC = 0f;
        public int PkCount = 0;
        public float ScanStartTime = 0f;
        public int MassResolvingPower = 0;
        public MS1Scan Next = null;

        public void Splat()
        {
            MS1Scan Runner = this.Next;
            while (Runner != null)
            {
                Console.WriteLine(Runner.TIC + "\t" + Runner.PkCount + "\t" + Runner.ScanStartTime + "\t" + Runner.MassResolvingPower);
                Runner = Runner.Next;
            }
        }
    }

    class ScanMetrics
    {
        public string NativeID = "";
        public int ScanNumber = 0;
        public float ScanStartTime = 0f;
        public int mzMLPeakCount = 0;
        public int MassResolvingPower = 0;
        public float MobilityValue = 0f;
        public float mzMLtic;
        public double IsolationTarget;
        public double IsolationLowerOffset;
        public double IsolationHigherOffset;
        public ScanMetrics Next;
        //TODO Allow users to specify alternative mass for Cys
        // Values from https://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html
        public static double[] AminoAcids = {57.02146,71.03711,87.03203,97.05276,99.06841,101.04768,103.00919,113.08406,114.04293,
                          115.02694,128.05858,128.09496,129.04259,131.04049,137.05891,147.06841,156.10111,
                          163.06333,186.07931};
        public static string[] AminoAcidSymbols = { "G", "A", "S", "P", "V", "T", "C", "L/I", "N",
                            "D", "Q", "K", "E", "M", "H", "F", "R", "Y", "W" };

        public ScanMetrics ExciseAllAtFirstIsolationTarget()
        {
            // Separate all scans that match the first isolation window in the current list.  Return that sublist.
            ScanMetrics SMRunner = this.Next;
            if (SMRunner == null) return null;
            double TargetMZ = SMRunner.IsolationTarget;
            float TargetFAIMS = SMRunner.MobilityValue;
            // Console.WriteLine("Extracting an isolation window at " + TargetMZ.ToString() + " with FAIMS CV of " + TargetFAIMS.ToString());
            int HowMany = 0;
            ScanMetrics SMHeader = new ScanMetrics();
            ScanMetrics SMRunner2 = SMHeader;
            SMRunner = this;
            while (SMRunner != null)
            {
                if (SMRunner.Next != null)
                {
                    if ((SMRunner.Next.MobilityValue == TargetFAIMS) && (SMRunner.Next.IsolationTarget == TargetMZ))
                    {
                        HowMany++;
                        SMRunner2.Next = SMRunner.Next;
                        SMRunner2 = SMRunner2.Next;
                        SMRunner.Next = SMRunner2.Next;
                        SMRunner2.Next = null;
                    }
                    else SMRunner = SMRunner.Next;
                }
                else SMRunner = SMRunner.Next;
            }
            // Console.WriteLine("\t" + HowMany.ToString() + " MS/MS measurements collected.");
            return SMHeader;
        }

        public string Stringify()
        {
            return (this.MobilityValue.ToString() + " " + this.IsolationTarget.ToString() + " " + this.ScanStartTime.ToString());
        }

        public bool ComesBefore(ScanMetrics Other)
        {
            if (this.MobilityValue < Other.MobilityValue)
            {
                return true;
            }
            if (this.IsolationTarget < Other.IsolationTarget)
            {
                return true;
            }
            if (this.ScanStartTime < Other.ScanStartTime)
            {
                return true;
            }
            return false;
        }
    }

    class SWATHMetrics
    {
        public double LoMZ = 0;
        public double HiMZ = 0;
        public double TargetMZ = 0;
        public double WidthMZ = 0;
        public float FAIMS = 0;
        public int MassResolvingPower = 0;
        public int MSMSCount = 0;
        public float LoRT = 0;
        public float HiRT = 0;
        public float CycleTimeMedian = 0;
        public float TIC25ileRT = 0;
        public float TIC50ileRT = 0;
        public float TIC75ileRT = 0;
        public float TotalTIC = 0;
        public int PkCountMin = 0;
        public int PkCount25ile = 0;
        public int PkCount50ile = 0;
        public int PkCount75ile = 0;
        public int PkCountMax = 0;
        public SWATHMetrics Next;
    }

    class LCMSMSExperiment
    {
        // Fields read directly from file
        public string SourceFile = "";
        public Uri FileURI;
        public string Instrument = "";
        public string SerialNumber = "";
        public string StartTimeStamp = "";
        public float MaxScanStartTime;
        // Computed fields
        public int mzMLMS1Count;
        public int mzMLMSnCount;
        public int MS1Resolution = 0;
        public float MS1TIC25ileRT = 0f;
        public float MS1TIC50ileRT = 0f;
        public float MS1TIC75ileRT = 0f;
        public float MS1TotalTIC = 0f;
        public float MS1CycleTime = 0f;
        public int MS1PkCountMin = 0;
        public int MS1PkCount25ile = 0;
        public int MS1PkCount50ile = 0;
        public int MS1PkCount75ile = 0;
        public int MS1PkCountMax = 0;
        public int SWATHCount = 0;
        public int SWATHCycleCountMin = Int32.MaxValue;
        public int SWATHCycleCountMax = 0;
        public double LoMZRange = Double.PositiveInfinity;
        public double HiMZRange = 0f;
        public double SWATHWidest = 0f;
        public double SWATHNarrowest = Double.PositiveInfinity;
        public float MedianCycleTimeAverage = 0f;
        public float TIC50ileRTMin = Single.PositiveInfinity;
        public float TIC50ileRTMax = 0f;
        public float TotalTICMin = Single.PositiveInfinity;
        public float TotalTICMax = 0f;
        public int PkCount50ileMin = Int32.MaxValue;
        public int PkCount50ileMax = 0;
        public static int MaxPkCount = 10000;
        public int[] mzMLPeakCountDistn = new int[MaxPkCount + 1];
        // Per-scan metrics
        public MS1Scan MS1Table = new MS1Scan();
        public ScanMetrics ScansTable = new ScanMetrics();
        public SWATHMetrics SWATHs = new SWATHMetrics();
        private MS1Scan MS1Runner;
        private ScanMetrics ScansRunner;
        public LCMSMSExperiment Next;
        private int LastPeakCount;
        private string LastNativeID = "";

        public void ComputeMS1Metrics()
        {
            MS1Runner = this.MS1Table.Next;
            float TICSum = 0f;
            float[] TICArray = new float[this.mzMLMS1Count];
            float[] RTArray = new float[this.mzMLMS1Count];
            int[] PkCountArray = new int[this.mzMLMS1Count];
            int Counter = 0;
            while (MS1Runner != null)
            {
                TICSum += MS1Runner.TIC;
                TICArray[Counter] = MS1Runner.TIC;
                RTArray[Counter] = MS1Runner.ScanStartTime;
                PkCountArray[Counter] = MS1Runner.PkCount;
                MS1Resolution = MS1Runner.MassResolvingPower;
                Counter++;
                MS1Runner = MS1Runner.Next;
            }
            this.MS1TotalTIC = TICSum;
            float TIC25 = TICSum / 4;
            float TIC50 = TICSum / 2;
            float TIC75 = TIC25 + TIC50;
            float TICSoFar = 0;
            float TICAfterThisScan = 0;
            int PkCountIndex = 0;
            if (this.mzMLMS1Count > 0)
            {
                float[] CycleTimes = new float[this.mzMLMS1Count - 1];
                float LastScanStartTime = 0;
                MS1Runner = this.MS1Table.Next;
                while (MS1Runner != null)
                {
                    TICAfterThisScan = TICSoFar + MS1Runner.TIC;
                    if ((TIC25 > TICSoFar) && (TIC25 <= TICAfterThisScan)) this.MS1TIC25ileRT = MS1Runner.ScanStartTime;
                    if ((TIC50 > TICSoFar) && (TIC50 <= TICAfterThisScan)) this.MS1TIC50ileRT = MS1Runner.ScanStartTime;
                    if ((TIC75 > TICSoFar) && (TIC75 <= TICAfterThisScan)) this.MS1TIC75ileRT = MS1Runner.ScanStartTime;
                    TICSoFar = TICAfterThisScan;
                    if (PkCountIndex > 0) CycleTimes[PkCountIndex - 1] = MS1Runner.ScanStartTime - LastScanStartTime;
                    PkCountIndex++;
                    LastScanStartTime = MS1Runner.ScanStartTime;
                    MS1Runner = MS1Runner.Next;
                }
                Array.Sort(TICArray);
                Array.Sort(PkCountArray);
                Array.Sort(CycleTimes);
                this.MS1PkCountMin = PkCountArray[0];
                this.MS1PkCount25ile = PkCountArray[this.mzMLMS1Count / 4];
                this.MS1PkCount50ile = PkCountArray[this.mzMLMS1Count / 2];
                this.MS1PkCount75ile = PkCountArray[(this.mzMLMS1Count / 4) + (this.mzMLMS1Count / 2)];
                this.MS1PkCountMax = PkCountArray[this.mzMLMS1Count - 1];
                //We multiply by 60 to get seconds from minutes.
                this.MS1CycleTime = 60 * CycleTimes[this.mzMLMS1Count / 2];
            }
        }

        public void ComputeMetricsForSwaths()
        {
            ScanMetrics TheseScans;
            ScanMetrics SMRunner;
            SWATHMetrics SWATHRunner = this.SWATHs;
            TheseScans = this.ScansTable.ExciseAllAtFirstIsolationTarget();
            float SumMedianCycleTime = 0;
            while (TheseScans != null)
            {
                //Compute SWATH metrics
                SWATHRunner.Next = new SWATHMetrics();
                SWATHRunner = SWATHRunner.Next;
                this.SWATHCount++;
                // We live dangerously, assuming ScanMetrics are ordered in retention time
                SMRunner = TheseScans.Next;
                SWATHRunner.LoRT = SMRunner.ScanStartTime;
                SWATHRunner.TargetMZ = SMRunner.IsolationTarget;
                SWATHRunner.LoMZ = SMRunner.IsolationTarget - SMRunner.IsolationLowerOffset;
                SWATHRunner.HiMZ = SMRunner.IsolationTarget + SMRunner.IsolationHigherOffset;
                SWATHRunner.WidthMZ = SWATHRunner.HiMZ - SWATHRunner.LoMZ;
                SWATHRunner.FAIMS = SMRunner.MobilityValue;
                SWATHRunner.MassResolvingPower = SMRunner.MassResolvingPower;
                if (SWATHRunner.LoMZ < this.LoMZRange) this.LoMZRange = SWATHRunner.LoMZ;
                if (SWATHRunner.HiMZ > this.HiMZRange) this.HiMZRange = SWATHRunner.HiMZ;
                if (SWATHRunner.WidthMZ > this.SWATHWidest) this.SWATHWidest = SWATHRunner.WidthMZ;
                if (SWATHRunner.WidthMZ < this.SWATHNarrowest) this.SWATHNarrowest = SWATHRunner.WidthMZ;
                float TICSum = 0;
                while (SMRunner != null)
                {
                    SWATHRunner.MSMSCount++;
                    SWATHRunner.HiRT = SMRunner.ScanStartTime;
                    TICSum = TICSum + SMRunner.mzMLtic;
                    if (SMRunner.mzMLPeakCount > SWATHRunner.PkCountMax) SWATHRunner.PkCountMax = SMRunner.mzMLPeakCount;
                    SMRunner = SMRunner.Next;
                }
                SWATHRunner.TotalTIC = TICSum;
                if (SWATHRunner.TotalTIC < this.TotalTICMin) this.TotalTICMin = SWATHRunner.TotalTIC;
                if (SWATHRunner.TotalTIC > this.TotalTICMax) this.TotalTICMax = SWATHRunner.TotalTIC;
                SMRunner = TheseScans.Next;
                float TIC25 = TICSum / 4;
                float TIC50 = TICSum / 2;
                float TIC75 = TIC25 + TIC50;
                float TICSoFar = 0;
                float TICAfterThisScan = 0;
                int[] PkCounts = new int[SWATHRunner.MSMSCount];
                int PkCountIndex = 0;
                float[] CycleTimes = new float[SWATHRunner.MSMSCount - 1];
                float LastScanStartTime = 0;
                if (SWATHRunner.MSMSCount < this.SWATHCycleCountMin) this.SWATHCycleCountMin = SWATHRunner.MSMSCount;
                if (SWATHRunner.MSMSCount > this.SWATHCycleCountMax) this.SWATHCycleCountMax = SWATHRunner.MSMSCount;
                while (SMRunner != null)
                {
                    TICAfterThisScan = TICSoFar + SMRunner.mzMLtic;
                    if ((TIC25 > TICSoFar) && (TIC25 <= TICAfterThisScan)) SWATHRunner.TIC25ileRT = SMRunner.ScanStartTime;
                    if ((TIC50 > TICSoFar) && (TIC50 <= TICAfterThisScan)) SWATHRunner.TIC50ileRT = SMRunner.ScanStartTime;
                    if ((TIC75 > TICSoFar) && (TIC75 <= TICAfterThisScan)) SWATHRunner.TIC75ileRT = SMRunner.ScanStartTime;
                    TICSoFar = TICAfterThisScan;
                    PkCounts[PkCountIndex] = SMRunner.mzMLPeakCount;
                    if (PkCountIndex > 0) CycleTimes[PkCountIndex - 1] = SMRunner.ScanStartTime - LastScanStartTime;
                    PkCountIndex++;
                    LastScanStartTime = SMRunner.ScanStartTime;
                    SMRunner = SMRunner.Next;
                }
                if (SWATHRunner.TIC50ileRT < this.TIC50ileRTMin) this.TIC50ileRTMin = SWATHRunner.TIC50ileRT;
                if (SWATHRunner.TIC50ileRT > this.TIC50ileRTMax) this.TIC50ileRTMax = SWATHRunner.TIC50ileRT;
                Array.Sort(PkCounts);
                SWATHRunner.PkCountMin = PkCounts[0];
                SWATHRunner.PkCount25ile = PkCounts[SWATHRunner.MSMSCount / 4];
                SWATHRunner.PkCount50ile = PkCounts[SWATHRunner.MSMSCount / 2];
                SWATHRunner.PkCount75ile = PkCounts[(SWATHRunner.MSMSCount / 4) + (SWATHRunner.MSMSCount / 2)];
                SWATHRunner.PkCountMax = PkCounts[SWATHRunner.MSMSCount - 1];
                if (SWATHRunner.PkCount50ile < this.PkCount50ileMin) this.PkCount50ileMin = SWATHRunner.PkCount50ile;
                if (SWATHRunner.PkCount50ile > this.PkCount50ileMax) this.PkCount50ileMax = SWATHRunner.PkCount50ile;
                Array.Sort(CycleTimes);
                //We multiply by 60 to get seconds from minutes.
                SWATHRunner.CycleTimeMedian = 60 * CycleTimes[SWATHRunner.MSMSCount / 2];
                SumMedianCycleTime += SWATHRunner.CycleTimeMedian;
                TheseScans = ScansTable.ExciseAllAtFirstIsolationTarget();
            }
            this.MedianCycleTimeAverage = SumMedianCycleTime / (float)this.SWATHCount;
        }

        public void BubbleSortScans()
        {
            /* Our goal is to get all the MS/MS scans of the same type
             * sorted together, ordered by their retention times.
             * Ideally, they are already sorted by retention time, but
             * it's possible that they are not. */
            ScanMetrics SMRunner;
            ScanMetrics SMBuffer1;
            ScanMetrics SMBuffer2;
            bool MadeChanges = true;
            while (MadeChanges)
            {
                MadeChanges = false;
                SMRunner = this.ScansTable;
                while (SMRunner != null)
                {
                    SMBuffer1 = SMRunner.Next;
                    if (SMBuffer1 != null)
                    {
                        SMBuffer2 = SMBuffer1.Next;
                        if (SMBuffer2 != null)
                        {
                            if (SMBuffer2.ComesBefore(SMBuffer1))
                            {
                                Console.WriteLine("Swapping these two:");
                                Console.WriteLine(SMBuffer1.Stringify());
                                Console.WriteLine(SMBuffer2.Stringify());
                                // Reorder the pair after SMRunner
                                SMRunner.Next = SMBuffer2;
                                SMBuffer1.Next = SMBuffer2.Next;
                                SMBuffer2.Next = SMBuffer1;
                                MadeChanges = true;
                            }
                        }
                    }
                    SMRunner = SMRunner.Next;
                }
            }
        }

        public void ReadFromMZML(XmlReader Xread)
        {
            /*
              Do not think of this code as a general-purpose mzML
              reader.  It is intended to populate only the fields that
              DIAuditor cares about.  For example, it entirely ignores
              the array of m/z values and intensities stored for any
              spectra.  This is intended to glean only the required
              fields in a single pass of the file.  It uses only the
              System.Xml libraries from Microsoft, obviating the need
              for any add-in libraries (to simplify the build process
              to something even I can use).
             */
            ScansRunner = ScansTable;
            MS1Runner = MS1Table;
            int CurrentScanType = 1;
            while (Xread.Read())
            {
                var ThisNodeType = Xread.NodeType;
                if (ThisNodeType == XmlNodeType.Element)
                {
                    if (Xread.Name == "run")
                    {
                        /*
                          We directly read relatively little
                          information about the mzML file as a whole.
                          Here we grab the startTimeStamp, and
                          elsewhere we grab the file name root,
                          instrument model, and serial number.
                        */
                        StartTimeStamp = Xread.GetAttribute("startTimeStamp");
                    }
                    else if (Xread.Name == "spectrum")
                    {
                        /*
                          We only create a new ScanMetrics object if
                          it isn't an MS1 scan.  We need to keep two
                          pieces of information from this new spectrum
                          header in case we do make a new ScanMetrics
                          object.
                        */
                        var ThisPeakCount = Xread.GetAttribute("defaultArrayLength");
                        LastPeakCount = int.Parse(ThisPeakCount);
                        LastNativeID = Xread.GetAttribute("id");
                    }
                    else if (Xread.Name == "cvParam")
                    {
                        var Accession = Xread.GetAttribute("accession");
                        switch (Accession)
                        {
                            /*
                              If you see that instrument model is ever
                              blank for an mzML, there are two likely
                              causes.  The first would be that the
                              mzML converter has not listed the CV
                              term for the instrument type in the
                              mzML-- ProteoWizard _does_ record this
                              information.  The second likely cause is
                              that the CV term relating to your
                              instrument model is missing from this
                              list.  Just add a "case" line for it and
                              recompile.
                             */
                            case "MS:1000126":
                            case "MS:1000557":
                            case "MS:1000932":
                            case "MS:1001742":
                            case "MS:1001910":
                            case "MS:1001911":
                            case "MS:1002416":
                            case "MS:1002523":
                            case "MS:1002533":
                            case "MS:1002634":
                            case "MS:1002732":
                            case "MS:1002877":
                            case "MS:1003005":
                            case "MS:1003028":
                            case "MS:1003029":
                            case "MS:1003094":
                            case "MS:1003096":
                            case "MS:1003123":
                            case "MS:1003293":
                            case "MS:1003356":
                            case "MS:1003378":
                                Instrument = Xread.GetAttribute("name");
                                break;
                            case "MS:1000529":
                                SerialNumber = Xread.GetAttribute("value");
                                break;
                            case "MS:1000016":
                                var ThisStartTime = Xread.GetAttribute("value");
                                // We need the "InvariantCulture" nonsense because some parts of the world separate decimals with commas.
                                var ThisStartTimeFloat = Single.Parse(ThisStartTime, CultureInfo.InvariantCulture);
                                if (Xread.GetAttribute("unitAccession")== "UO:0000031")
                                {
                                    // We are in minutes; do nothing
                                }
                                else if (Xread.GetAttribute("unitAccession")== "UO:0000028")
                                {
                                    // We are in seconds; divide by 60.
                                    ThisStartTimeFloat = ThisStartTimeFloat / 60.0f;
                                }
                                else
                                {
                                    WarningException myEx = new WarningException("scan start times are written in something other than minutes or seconds!");
                                    Console.Write(myEx.ToString());
                                }
                                if (CurrentScanType == 1)
                                {
                                    MS1Runner.ScanStartTime = ThisStartTimeFloat;
                                }
                                else
                                {
                                    ScansRunner.ScanStartTime = ThisStartTimeFloat;
                                    if (ThisStartTimeFloat > MaxScanStartTime) MaxScanStartTime = ThisStartTimeFloat;
                                }
                                break;
                            case "MS:1001581":
                                var ThisFAIMS = Xread.GetAttribute("value");
                                ScansRunner.MobilityValue = Single.Parse(ThisFAIMS, CultureInfo.InvariantCulture);
                                break;
                            case "MS:1000827":
                                var ThisIsolationTarget = Xread.GetAttribute("value");
                                ScansRunner.IsolationTarget = Double.Parse(ThisIsolationTarget, CultureInfo.InvariantCulture);
                                break;
                            case "MS:1000828":
                                var ThisIsolationLower = Xread.GetAttribute("value");
                                ScansRunner.IsolationLowerOffset = Double.Parse(ThisIsolationLower, CultureInfo.InvariantCulture);
                                break;
                            case "MS:1000829":
                                var ThisIsolationHigher = Xread.GetAttribute("value");
                                ScansRunner.IsolationHigherOffset = Double.Parse(ThisIsolationHigher, CultureInfo.InvariantCulture);
                                break;
                            case "MS:1000285":
                                var ThisTIC = Xread.GetAttribute("value");
                                if (CurrentScanType == 1)
                                {
                                    MS1Runner.TIC = Single.Parse(ThisTIC, CultureInfo.InvariantCulture);
                                }
                                else
                                {
                                    ScansRunner.mzMLtic = Single.Parse(ThisTIC, CultureInfo.InvariantCulture);
                                }
                                break;
                            case "MS:1000800":
                                var ThisPower = Xread.GetAttribute("value");
                                if (CurrentScanType == 1)
                                {
                                    MS1Runner.MassResolvingPower = int.Parse(ThisPower);
                                }
                                else
                                {
                                    ScansRunner.MassResolvingPower = int.Parse(ThisPower);
                                }
                                break;
                            case "MS:1000511":
                                var ThisLevel = Xread.GetAttribute("value");
                                CurrentScanType = int.Parse(ThisLevel);
                                if (CurrentScanType == 1)
                                {
                                    mzMLMS1Count++;
                                    MS1Runner.Next = new MS1Scan();
                                    MS1Runner = MS1Runner.Next;
                                    MS1Runner.PkCount = LastPeakCount;
                                }
                                else
                                {
                                    /*
                                      If we detect an MS of level 2 or
                                      greater in the mzML, we have
                                      work to do.  Each MS/MS is
                                      matched by an item in the linked
                                      list of ScanMetrics for each
                                      LCMSMSExperiment.  We will need
                                      to capture some information we
                                      already saw (such as the
                                      NativeID of this scan) and set
                                      up for collection of some
                                      additional information in the
                                      Activation section.
                                     */
                                    mzMLMSnCount++;
                                    ScansRunner.Next = new ScanMetrics();
                                    ScansRunner = ScansRunner.Next;
                                    ScansRunner.NativeID = LastNativeID;
                                    ScansRunner.mzMLPeakCount = LastPeakCount;
                                    if (LastPeakCount > MaxPkCount)
                                        mzMLPeakCountDistn[MaxPkCount]++;
                                    else
                                        mzMLPeakCountDistn[LastPeakCount]++;
                                }
                                break;
                        }
                    }
                }
            }
        }

        public LCMSMSExperiment Find(string Basename)
        {
            /*
              We receive an msAlign filename.  We seek the
              LCMSMSExperiment in this linked list that has the
              corresponding filename.
            */
            var LRunner = this.Next;
            while (LRunner != null)
            {
                if (Basename == LRunner.SourceFile)
                    return LRunner;
                LRunner = LRunner.Next;
            }
            return null;
        }

        public void ParseScanNumbers()
        {
            /*
              Instrument manufacturers differ in the ways that they
              report the identities of each MS and MS/MS.  Because
              TopFD was created for Thermo instruments, though, it
              expects each spectrum to have a unique scan number for a
              given RAW file.  To match to TopFD msAlign files, we'll
              need to extract those from the NativeIDs.
             */
            var SRunner = this.ScansTable.Next;
            while (SRunner != null)
            {
                // Example of Thermo NativeID: controllerType=0 controllerNumber=1 scan=12 (ProteoWizard)
                // Example of SCIEX NativeID: sample=1 period=1 cycle=806 experiment=2 (ProteoWizard)
                // Example of SCIEX NativeID: sample=1 period=1 cycle=7207 experiment=1 (SCIEX MS Data Converter)
                // Example of Bruker NativeID: scan=55 (TIMSConvert)
                // Example of Bruker NativeID: merged=102 frame=13 scanStart=810 scanEnd=834 (ProteoWizard)
                var Tokens = SRunner.NativeID.Split(' ');
                foreach (var ThisTerm in Tokens)
                {
                    var Tokens2 = ThisTerm.Split('=');
                    if (Tokens2[0].Equals("cycle") || Tokens2[0].Equals("scan") || Tokens2[0].Equals("scanStart"))
                    {
                        SRunner.ScanNumber = int.Parse(Tokens2[1]);
                    }
                }
                SRunner = SRunner.Next;
            }
        }

        public ScanMetrics GoToScan(int Target)
        {
            var SRunner = this.ScansTable.Next;
            while (SRunner != null)
            {
                if (SRunner.ScanNumber == Target)
                    return SRunner;
                SRunner = SRunner.Next;
            }
            return null;
        }

        public int[] QuartilesOf(int[] Histogram)
        {
            var Quartiles = new int[5];
            var Sum = 0;
            int index;
            var AwaitingMin = true;
            for (index = 0; index < Histogram.Length; index++)
            {
                if (AwaitingMin && Histogram[index] > 0)
                {
                    AwaitingMin = false;
                    Quartiles[0] = index;
                }
                if (Histogram[index] > 0)
                    Quartiles[4] = index;
                Sum += Histogram[index];
            }
            var CountQ1 = Sum / 4;
            var CountQ2 = Sum / 2;
            var CountQ3 = CountQ1 + CountQ2;
            Sum = 0;
            for (index = 0; index < Histogram.Length; index++)
            {
                var ThisCount = Histogram[index];
                if (Sum < CountQ1 && CountQ1 <= Sum + ThisCount)
                    Quartiles[1] = index;
                if (Sum < CountQ2 && CountQ2 <= Sum + ThisCount)
                    Quartiles[2] = index;
                if (Sum < CountQ3 && CountQ3 <= Sum + ThisCount)
                    Quartiles[3] = index;
                Sum += ThisCount;
            }
            return Quartiles;
        }
        public void WriteMZQCReport(string Version)
        {
            var LCMSMSRunner = this.Next;
            var cv = new ControlledVocabulary
            {
                Name = "Proteomics Standards Initiative Mass Spectrometry Ontology",
                Uri = new System.Uri("https://github.com/HUPO-PSI/psi-ms-CV/releases/download/v4.1.186/psi-ms.obo"),
                Version = "4.1.186"
            };
            var runs = new List<BaseQuality>();
            while (LCMSMSRunner != null)
            {
                float[] MS1TICQuartileRTs = { LCMSMSRunner.MS1TIC25ileRT, LCMSMSRunner.MS1TIC50ileRT, LCMSMSRunner.MS1TIC75ileRT };
                double[] SWATHmzRange = { LCMSMSRunner.LoMZRange, LCMSMSRunner.HiMZRange };
                double[] SWATHWidthRange = { LCMSMSRunner.SWATHNarrowest, LCMSMSRunner.SWATHWidest };
                float[] SWATHRetentionTimeRangeAtHalfTIC = { LCMSMSRunner.TIC50ileRTMin, LCMSMSRunner.TIC50ileRTMax };
                float[] SWATHTICSumRange = { LCMSMSRunner.TotalTICMin, LCMSMSRunner.TotalTICMax };
                int[] SWATHPkCountMedianRange = { LCMSMSRunner.PkCount50ileMin, LCMSMSRunner.PkCount50ileMax };
                int[] MS1PkCountQuartiles = { LCMSMSRunner.MS1PkCountMin, LCMSMSRunner.MS1PkCount25ile, LCMSMSRunner.MS1PkCount50ile, LCMSMSRunner.MS1PkCount75ile, LCMSMSRunner.MS1PkCountMax };
                int[] SWATHCycleCountRange = { LCMSMSRunner.SWATHCycleCountMin, LCMSMSRunner.SWATHCycleCountMax };
                //Why does FileFormat write a null "value" field?
                var run = new BaseQuality
                {
                    Metadata = new Metadata
                    {
                        Label = LCMSMSRunner.SourceFile,
                        AnalysisSoftware = [new AnalysisSoftwareElement { Accession = "MS:4000189", Name = "DIAuditor",
                            Value = "Data-Independent Acquisition QC Metric Generator", Uri = new Uri("https://github.com/dtabb73/DIAuditor"), Version=Version}],
                        InputFiles = [ new InputFile { Name = LCMSMSRunner.SourceFile,
                            Location = LCMSMSRunner.FileURI,
                            FileFormat = new CvParameter { Accession = "MS:1000584", Name = "mzML format"},
                            FileProperties = [new CvParameter { Accession = "MS:1000031",
                                Name = "instrument model",
                                Value = LCMSMSRunner.Instrument},
                                new CvParameter { Accession = "MS:1000747",
                                Name = "completion time",
                                Value = LCMSMSRunner.StartTimeStamp} ]
                        }]
                    }
                };
                run.QualityMetrics = new List<QualityMetric>();
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:1000529",
                    Name = "instrument serial number",
                    Value = LCMSMSRunner.SerialNumber
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000053",
                    Name = "chromatography duration",
                    Value = LCMSMSRunner.MaxScanStartTime
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000059",
                    Name = "number of MS1 spectra",
                    Value = LCMSMSRunner.mzMLMS1Count
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000060",
                    Name = "number of MS2 spectra",
                    Value = LCMSMSRunner.mzMLMSnCount
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    //Can I use the same CV term for MS1 resolution as I use for MS2 resolution?
                    Accession = "MS:1000800",
                    Name = "MS1 mass resolving power",
                    Value = LCMSMSRunner.MS1Resolution
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000190",
                    Name = "The retention times at which 25%, 50%, and 75% of MS1 TIC have been accumulated",
                    Value = MS1TICQuartileRTs
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:1000285",
                    Name = "Summed Total Ion Current from MS1",
                    Value = LCMSMSRunner.MS1TotalTIC
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000192",
                    Name = "The median cycle time between MS1 measurements",
                    Value = LCMSMSRunner.MS1CycleTime
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000061",
                    Name = "MS1 peak count quartiles",
                    Value = MS1PkCountQuartiles
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000194",
                    Name = "Count of isolation windows in this DIA method",
                    Value = LCMSMSRunner.SWATHCount
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000196",
                    Name = "Range of MS/MS measurements per isolation window in this DIA method",
                    Value = SWATHCycleCountRange
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:1003159",
                    Name = "Range of lowest m/z in lowest isolation window to highest m/z of highest isolation window in this DIA method",
                    Value = SWATHmzRange
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000195",
                    Name = "Range of isolation window widths in this DIA method",
                    Value = SWATHWidthRange
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000197",
                    Name = "The range of retention times at half TIC among isolation windows",
                    Value = SWATHRetentionTimeRangeAtHalfTIC
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000198",
                    Name = "The range of TIC sums among isolation windows",
                    Value = SWATHTICSumRange
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000199",
                    Name = "The range of Peak Count medians among isolation windows",
                    Value = SWATHPkCountMedianRange
                }
                );
                // Now add metrics corresponding to each isolation window separately.
                double[] IsolationWindowTargets = new double[LCMSMSRunner.SWATHCount];
                double[] IsolationWindowOffsetLo = new double[LCMSMSRunner.SWATHCount];
                double[] IsolationWindowOffsetHi = new double[LCMSMSRunner.SWATHCount];
                float[] MobilityValues = new float[LCMSMSRunner.SWATHCount];
                int[] ResolvingPowers = new int[LCMSMSRunner.SWATHCount];
                int[] MSMSCounts = new int[LCMSMSRunner.SWATHCount];
                float[] MedianCycleTimes = new float[LCMSMSRunner.SWATHCount];
                float[,] RetentionTimesAtQuartilesOfTIC = new float[LCMSMSRunner.SWATHCount,3];
                float[] TotalTICs = new float[LCMSMSRunner.SWATHCount];
                int[,] PeakCountQuartiles = new int[LCMSMSRunner.SWATHCount,5];
                SWATHMetrics SMRunner = LCMSMSRunner.SWATHs.Next;
                int index = 0;
                while (SMRunner != null)
                {
                    IsolationWindowTargets[index] = SMRunner.TargetMZ;
                    IsolationWindowOffsetLo[index] = SMRunner.TargetMZ - SMRunner.LoMZ;
                    IsolationWindowOffsetHi[index] = SMRunner.HiMZ - SMRunner.TargetMZ;
                    MobilityValues[index] = SMRunner.FAIMS;
                    ResolvingPowers[index] = SMRunner.MassResolvingPower;
                    MSMSCounts[index] = SMRunner.MSMSCount;
                    MedianCycleTimes[index] = SMRunner.CycleTimeMedian;
                    RetentionTimesAtQuartilesOfTIC[index, 0] = SMRunner.TIC25ileRT;
                    RetentionTimesAtQuartilesOfTIC[index, 1] = SMRunner.TIC50ileRT;
                    RetentionTimesAtQuartilesOfTIC[index, 2] = SMRunner.TIC75ileRT;
                    TotalTICs[index] = SMRunner.TotalTIC;
                    PeakCountQuartiles[index, 0] = SMRunner.PkCountMin;
                    PeakCountQuartiles[index, 1] = SMRunner.PkCount25ile;
                    PeakCountQuartiles[index, 2] = SMRunner.PkCount50ile;
                    PeakCountQuartiles[index, 3] = SMRunner.PkCount75ile;
                    PeakCountQuartiles[index, 4] = SMRunner.PkCountMax;
                    index++;
                    SMRunner = SMRunner.Next;
                }
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:1000827",
                    Name = "The target m/z values for each isolation window",
                    Value = IsolationWindowTargets
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:1000828",
                    Name = "The lower offset m/z values for each isolation window",
                    Value = IsolationWindowOffsetLo
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:1000829",
                    Name = "The upper offset m/z values for each isolation window",
                    Value = IsolationWindowOffsetHi
                }
                );
                run.QualityMetrics.Add(new QualityMetric
                {
                    // We want this to provide information on tims thresholds, too, but that's not currently included in mzML data model.
                    Accession = "MS:1001581",
                    Name = "FAIMS compensation voltage for each isolation window",
                    Value = MobilityValues
                });
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:1000234",
                    Name = "Mass Resolving Power for each isolation window",
                    Value = ResolvingPowers
                });
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000060",
                    Name = "Number of MS2 Spectra for each isolation window",
                    Value = MSMSCounts
                });
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000193",
                    Name = "The median cycle time between MS2 measurements for each isolation window",
                    Value = MedianCycleTimes
                });
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000191",
                    Name = "The retention times at which 25%, 50%, and 75% of MS2 TIC have been accumulated for each isolation window",
                    Value = RetentionTimesAtQuartilesOfTIC
                });
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:1000285",
                    Name = "Summed Total Ion Current from MS2 for each isolation window",
                    Value = TotalTICs
                });
                run.QualityMetrics.Add(new QualityMetric
                {
                    Accession = "MS:4000062",
                    Name = "The quartiles for peak counts per MS/MS for each isolation window",
                    Value = PeakCountQuartiles
                });
                runs.Add(run);
                LCMSMSRunner = LCMSMSRunner.Next;
            }
            var mzqc = new MzqcContent
            {
                Description = "DIAuditor QC report",
                Version = "1.0.0",
                CreationDate = DateTime.Now,
                ControlledVocabularies = [cv],
                RunQualities = runs
            };
            var file = new Mzqc { MzqcContent = mzqc };
            File.WriteAllText(@"DIAuditor.mzqc", file.ToJson());
        }

        public void WriteTextQCReport()
        {
            /*
              We have two TSV outputs.  The "byRun" report contains a
              row for each LC-MS/MS (or mzML) in this directory.  The
              "byMSn" report contains a row for each MS/MS in each
              mzML in this directory.
             */
            //TODO: Should I be reporting distribution of deconvolved precursor mass by RAW?
            var LCMSMSRunner = this.Next;
            const string delim = "\t";
            using (var TSVbyRun = new StreamWriter("DIAuditor-byRun.tsv"))
            {
                TSVbyRun.WriteLine("SourceFile\tInstrument\tSerialNumber\tStartTimeStamp\tRTDuration" +
                   "\tmzMLMS1Count\tmzMLMSnCount\tMS1Resolution" +
                   "\tMS1TIC25ileRT\tMS1TIC50ileRT\tMS1TIC75ileRT\tMS1TotalTIC\tMS1CycleTime" +
                   "\tMS1PkCountMin\tMS1PkCount25ile\tMS1PkCount50ile\tMS1PkCount75ile\tMS1PkCountMax" +
                   "\tIsolationWindowCount\tCyclesMin\tCyclesMax" +
                   "\tMZRangeMin\tMZRangeMax\tIsolationWindowWidthMin\tIsolationWidowWidthMax" +
                   "\tAverageMedianCycleTime\tTICMedianRTMin\tTICMedianRTMax\tTotalTICMin\tTotalTICMax" +
                   "\tPkCountMedianMin\tPkCountMedianMax");
                while (LCMSMSRunner != null)
                {
                    // Actually write the metrics to the byRun file...
                    TSVbyRun.Write(LCMSMSRunner.SourceFile + delim);
                    TSVbyRun.Write(LCMSMSRunner.Instrument + delim);
                    TSVbyRun.Write(LCMSMSRunner.SerialNumber + delim);
                    TSVbyRun.Write(LCMSMSRunner.StartTimeStamp + delim);
                    TSVbyRun.Write(LCMSMSRunner.MaxScanStartTime + delim);
                    TSVbyRun.Write(LCMSMSRunner.mzMLMS1Count + delim);
                    TSVbyRun.Write(LCMSMSRunner.mzMLMSnCount + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1Resolution + delim);

                    TSVbyRun.Write(LCMSMSRunner.MS1TIC25ileRT + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1TIC50ileRT + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1TIC75ileRT + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1TotalTIC + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1CycleTime + delim);

                    TSVbyRun.Write(LCMSMSRunner.MS1PkCountMin + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1PkCount25ile + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1PkCount50ile + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1PkCount75ile + delim);
                    TSVbyRun.Write(LCMSMSRunner.MS1PkCountMax + delim);

                    TSVbyRun.Write(LCMSMSRunner.SWATHCount + delim);
                    TSVbyRun.Write(LCMSMSRunner.SWATHCycleCountMin + delim);
                    TSVbyRun.Write(LCMSMSRunner.SWATHCycleCountMax + delim);
                    TSVbyRun.Write(LCMSMSRunner.LoMZRange + delim);
                    TSVbyRun.Write(LCMSMSRunner.HiMZRange + delim);
                    TSVbyRun.Write(LCMSMSRunner.SWATHNarrowest + delim);
                    TSVbyRun.Write(LCMSMSRunner.SWATHWidest + delim);
                    TSVbyRun.Write(LCMSMSRunner.MedianCycleTimeAverage + delim);
                    TSVbyRun.Write(LCMSMSRunner.TIC50ileRTMin + delim);
                    TSVbyRun.Write(LCMSMSRunner.TIC50ileRTMax + delim);
                    TSVbyRun.Write(LCMSMSRunner.TotalTICMin + delim);
                    TSVbyRun.Write(LCMSMSRunner.TotalTICMax + delim);
                    TSVbyRun.Write(LCMSMSRunner.PkCount50ileMin + delim);
                    TSVbyRun.WriteLine(LCMSMSRunner.PkCount50ileMax + delim);
                    /*
                            foreach (var ThisQuartile in LCMSMSRunner.mzMLPeakCountQuartiles)
                                TSVbyRun.Write(ThisQuartile + delim);
                    */
                    LCMSMSRunner = LCMSMSRunner.Next;
                }
            }
            LCMSMSRunner = this.Next;
            using (var TSVbySWATH = new StreamWriter("DIAuditor-byIsolationWindow.tsv"))
            {
                TSVbySWATH.WriteLine("SourceFile\tLoMZ\tHiMZ\tWidthMZ\tIonMobility\tMassResolvingPower" +
                    "\tMSMSCount\tRTMin\tRTMax\tCycleTimeMedian\tTIC25ileRT\tTIC50ileRT\tTIC75ileRT" +
                    "\tTotalTIC\tPkCountMin\tPkCount25ile\tPkCount50ile\tPkCount75ile\tPkCountMax");
                while (LCMSMSRunner != null)
                {
                    var SWATHRunner = LCMSMSRunner.SWATHs.Next;
                    while (SWATHRunner != null)
                    {
                        TSVbySWATH.Write(LCMSMSRunner.SourceFile + delim);
                        TSVbySWATH.Write(SWATHRunner.LoMZ + delim);
                        TSVbySWATH.Write(SWATHRunner.HiMZ + delim);
                        TSVbySWATH.Write(SWATHRunner.WidthMZ + delim);
                        TSVbySWATH.Write(SWATHRunner.FAIMS + delim);
                        TSVbySWATH.Write(SWATHRunner.MassResolvingPower + delim);
                        TSVbySWATH.Write(SWATHRunner.MSMSCount + delim);
                        TSVbySWATH.Write(SWATHRunner.LoRT + delim);
                        TSVbySWATH.Write(SWATHRunner.HiRT + delim);
                        TSVbySWATH.Write(SWATHRunner.CycleTimeMedian + delim);
                        TSVbySWATH.Write(SWATHRunner.TIC25ileRT + delim);
                        TSVbySWATH.Write(SWATHRunner.TIC50ileRT + delim);
                        TSVbySWATH.Write(SWATHRunner.TIC75ileRT + delim);
                        TSVbySWATH.Write(SWATHRunner.TotalTIC + delim);
                        TSVbySWATH.Write(SWATHRunner.PkCountMin + delim);
                        TSVbySWATH.Write(SWATHRunner.PkCount25ile + delim);
                        TSVbySWATH.Write(SWATHRunner.PkCount50ile + delim);
                        TSVbySWATH.Write(SWATHRunner.PkCount75ile + delim);
                        TSVbySWATH.WriteLine(SWATHRunner.PkCountMax);
                        SWATHRunner = SWATHRunner.Next;
                    }
                    LCMSMSRunner = LCMSMSRunner.Next;
                }
            }

        }
    }
}
