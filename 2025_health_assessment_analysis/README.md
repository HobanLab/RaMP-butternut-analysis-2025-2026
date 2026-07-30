# butternut-health-assessment-2025
This code analyzes the butternut health assessment (2025 version) data and seeks to understand the ecological conditions where butternuts survive and thrive.

General workflow I used:
1. Take data that has been cleaned in that it has been made correct and the associate trees are separated out.  The reason all of this is done in Sheets, not R, is because it mostly requires line-by-line changes that can't be caught by a computer (ex., unique typos, or data inputted incorrectly or in the wrong column) and judgement calls about how to group trees that also shouldn't be done by a computer.
2. Source DataSplitter, the code that creates smaller datasheets with subsets of data (ex. only the live adult trees, only the trees at a certain site) so that they can be analyzed individually.
3. Source ChartBuilder, the code that creates functions for commonly used graphics (ex. boxplots with significance testing, maps, etc.) so that they can be generated with one line of code, increasing efficiency and readability of code.
4. Run code that actually analyzes the data!  Answer these questions:
	a. Where are butternuts doing well?  Where are they regenerating? (mapping.R)
	b. How are butternut seedlings doing?  What factors contribute to successful seedling establishment? (seedlings.R)
	c. What ecological factors are associated with canker development and spread?
