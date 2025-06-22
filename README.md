# Organization
* renv files are now available to help load expected versions of dependencies. Use renv::install() or renv::restore() to get dependencies. Then follow Installation to install/build funclustweight.
* Files required for the package to run/build are inside of: R and src folders; NAMESPACE and DESCRIPTION files.
* Necessary documentation is inside of: inst (citations), LICENSE, and man (manual pages).
* Test files as well as necessary files for reproducability are inside of tests (a file named examples.R will run on build).
* Within the main directory a hidden file .RData contains an existing build with loaded packages and data.

# Beta Installation
1. Open RStudio
2. Select "File" > "New Project..."
3. In the dialog menu select "Existing Directory"
4. In project working directory select the path to funclustweight main directory and click "Create Project"
5. A new tab in the menu on the right should appear next to "Environment", "History", "Connections", etc. select "Build" from this menu
6. Click "Install" in the build menu and install any requirements (for a beta package one requirement will be RBuildTools, a menu will pop up to download these)
7. Run the line "library(funclustweight)" to load the built package, now you may run any required functions

# Help
* Run "?funclustweight" for documentation, there is also documentation for: fitAdelaideFD, genFromModelRegFD, plotAdelaideFD, and plotRegFD.
