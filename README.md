# LUX_C14_analysis
Author: Jun Liao (junhui.private@gmail.com, or junhui_liao@brown.edu).
This C-14 analysis aims to create a couple of acceptances of S1 signal for LUX Run4 data. The ER energy goes up to 350 phd which is required by EFT analysis. It is the extention of the previous CH3T analyiss which has an ER cutoff at 50 phd.
successCreateRootFiles.C is to access the full dataset file, "c14_endofLUX_09262017.root", to extract whatever interesting data, for instance, "S1_prompt_fraction_tlx" V.S. "S1_area", then save it as a ROOT file to be analyzed further.
analysisC14Data.C is to analyze "S1_prompt_fraction_tlx" V.S. "S1_area". In every bin of X-axis (S1_area), calculating the 99.5% bin content of its Y-axis (S1_prompt_fraction_tlx), to get an new array of [S1_are, 99.5% bin content of S1_prompt_fraction_tlx]. Fill this arry into a Graph and fitting with a combined function.
analyze_S1_max_peak_area_VS_S1_area.C is very similar to analysisC14Data.C. The difference is the Y-axis is "S1 max peak area".
