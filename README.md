# Global-Oil-Market
This code emulates the SVAR results from the paper: "Lutz Kilian, 2009. "Not All Oil Price Shocks Are Alike: Disentangling Demand and Supply Shocks in the Crude Oil Market," American Economic Review, vol. 99(3), pages 1053-1069, June.  Specifically:  Cholesky Decomposition, Structural Impulse Response, Historical Evolution of the Structural Shock,  Historical Decompositio of the Structural Shocks and the Forecast Error Variance Decomposition.



The sheet "Data_oil_1" contains monthly data for the period 1973:M1 to 2019:M6:

 RAC: US Refiner Acquisition Cost of Crude Oil (divided by US CPI) [in log]
 WTI: West Texas Intermediate Oil Price (divided by US CPI) [in log]
 Oil_Prod: Global Oil Production [in log]
 Oil_Inv: OECD Oil Inventories [in log]
 Kilian_Index: A proxy of global economic activity constructed from the freight rates by Lutz Kilian. [index]
 Hamilton_Index: 2 years growth in global industrial production [in %]
