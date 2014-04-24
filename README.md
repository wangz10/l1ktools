# L1000 Analysis Tools v1.0


Copyright 2011-2014 Broad Institute of MIT and Harvard.

A collection of software tools to read and analyze data produced from
the L1000 project (www.lincscloud.org).

## Analysis Tools

A brief description of the tools included in this software package is
given below. The Matlab implementation of the tools is currently the
most mature. Some basic utilities in R and java are also included. We
will update the tools as they become available.

### Matlab Tools: matlab/

#### Requirements:

1. Matlab R2009a and above
2. Statistics Toolbox

#### Tools:
* **dpeak_demo.m**: Demonstrates basic usage of the scripts. 
* **l1kt_dpeak.m**: Performs peak deconvolution for all analytes in a single LXB file, and outputs a report of the detected peaks.
* **l1kt_plot_peaks.m**: Plots intensity distributions for one or more analytes in an LXB file.
* **l1kt_parse_lxb.m**:	Reads an LXB file and returns the RID and RP1 values.
* **l1kt_liss.m**: Performs Luminex Invariant Set Smoothing on a raw (GEX) input .gct file
* **l1kt_qnorm.m**:	Performs quantile normalization on an input .gct file
* **l1kt_infer.m**:	Infers expression of target genes from expression of landmark genes in an input .gct file

See the documentation included with each script for a details on usage
and input parameters.

#### Demo:
* **dpeak_demo.m**: To run the demo, start Matlab, change to the folder containing dpeak_demo and
type dpeak_demo in the Command Window. This will read a sample LXB
file (A10.lxb), generate a number of intensity distribution plots and create a
text report of the statistics of the detected peaks (A10_pkstats.txt).

* **example_methods.m**: To run the demo, start Matlab, change to the folder containing example_methods and type example_methods at the command line. The script will read in a .gct and a .gctx file, z-score the data in the .gctx file, and read in an .lxb file.

### R Tools: R/

#### Requirements:

1. R versions 2.9 and above
2. prada package: http://www.bioconductor.org/packages/devel/bioc/html/prada.html
3. h5r package: http://cran.r-project.org/web/packages/h5r/index.html

#### Tools:

* **cmap/lxb2txt.R**	Saves values from an LXB file as a tab-delimited text file.
* **cmap/lxb2txt.sh** Bash wrapper to lxb2txt.R 

#### Demo:
* **example_methods.R**: To run the demo, change to the folder containing example_methods.R and source the script. It will read in a .gctx file and display its contents.


### Java Tools: java/

#### Tools:

* **lxb-util.jar**:	Java routine for parsing LXB files (source included).

### Data format

Bead-level measurements from each detected well in L1000 assay are
stored as binary .LXB files. The LXB files conform to the FCS 3.0 data
file standard, commonly used to store flow cytometry data.  For
details see: http://isac-net.org/home.aspx.

The two main measurement parameters are: 

* **RID**: 	<32 bit unsigned integer> Reporter identifier. Values ranging
	[1-500] specifies the bead color/region. RID=0 refer to events
	that flowcytometer software was unable to classify.

* **RP1**: 	<32 bit unsigned integer> Reporter fluorescent intensity,
	quantifies transcript abundance of the gene interrogated by
	the bead.
	
### Python Tools: python/

#### Requirements:

1. Python 2.7 (untested under Python 3)
2. numpy (http://numpy.scipy.org)
3. pandas (http://pandas.pydata.org/)
4. requests (http://docs.python-requests.org/en/latest/)
5. pytables (http://www.pytables.org/moin)
6. blessings (http://pypi.python.org/pypi/blessings)

#### Tools:
* **cmap/io/gct.py** : Classes to interact with .gct and .gctx files.
* **cmap/util/api_utils.py**: Classes to make calls to the LINCS annotation API and return results as Python data structures.

#### Demo:
* **example_methods.py**: To run the demo, change to the folder containing example_methods.py and run the script. It will read in a .gctx file and display its contents. It will then make calls to the LINCS web API to retrieve annotations. In order to make the API call, an API key must be provided. If you do not have a key, contact lincs@broadinstitute.org.

## Common data analysis tasks
Below are summarized the tools available to perform so common data analysis tasks.

### Reading .gct and .gctx files
* **MATLAB**: Add l1ktools/matlab to your path. Then use the parse_gctx function.
* **R**: Source the script l1ktools/R/cmap/io.R. Then use the parse.gctx function.
* **Python**: Add l1ktools/python to your Python path, and import the module cmap.io.gct. Then instantiate a GCT object, and call its read() method. For more information, see the documentation on the GCT class.

### Z-Scoring a data set
* **MATLAB**: Configure your MATLAB path as above. Then use the robust_zscore function.

### Interacting with the LINCS API
* **Python**: Configure your Python path as above, and import the module cmap.util.api_utils. The classes CMapAPI and APIContainer handle calls to the API; see their documentation for more details.

## Software License
For licensing information see [http://lincscloud.org/license/](http://lincscloud.org/license/).

