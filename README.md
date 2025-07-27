

[![DOI](https://zenodo.org/badge/612192520.svg)](https://zenodo.org/badge/latestdoi/612192520)


# GoldenStandardDataset
Workflow and datasets from the parts characterisation from Golden Standard Moclo kit.

In this repository you will find 3 jupyter notebooks:
- **paper_datasets_workflow**: this script allows you to explore our datasets and follow the same computational analysis we performed to quantify parts activity from growth and fluorescence data. You can find it here: <a target="_blank" href="https://colab.research.google.com/github/SBGlab/GoldenStandardDataset/blob/main/paper_datasets_workflow.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

- **custom_data_workflow**: script to analyse your own data following the same workflow. <a target="_blank" href="https://colab.research.google.com/github/SBGlab/GoldenStandardDataset/blob/main/custom_data_workflow.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a> You will need to previously generate two .csv files with a given structure*: as well as sample info file, which you can easily create using the sample_info_file script.

New v3: <a target="_blank" href="https://colab.research.google.com/github/SBGlab/GoldenStandardDataset/blob/main/create_sample_info_file_v3.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

- **sample_info_file**: script to easily generate sample-to-well info for processing custom experiments. We used this script also for each of our analysis. You can find it here: <a target="_blank" href="https://colab.research.google.com/github/SBGlab/GoldenStandardDataset/blob/main/create_sample_info_file.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

*The structure for .csv files is simple, leave first first cell empty and fill the first column with time data, first row with well data and each row with the measurements for that well and time point. Make sure you get a csv file separated by commas. Example files are provided within the custom example folder.

Note: You can upload the full 96 well plate data, unused wells will be automatically discarded if htey are not in the sample info file.

![alt text](https://github.com/SBGlab/GoldenStandardDataset/blob/main/csv_template.jpg)
