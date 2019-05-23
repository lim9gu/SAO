SAO1-m
=============
Data reduction package for Seoul National University Astronomical Observatory (SAO) 1-m telescope. Firstly, this code contains basic reduction of bias/dark subtraction, and flat fielding (skyflat, domeflat). Data classification (bias, Dark, domeflat, skyflat) can be performed. (Ex. Exposure time, Filter, Binning). Solving WCS is based on Astrometry.net software (Using 'Scamp' will be added soon). python2.7 and python3.6 is compatible.   

Current version is 0.5. This code will be updated continuously.   

Prerequisites
-------------  
**Python libraries :**  
astropy  
numpy  
Pyraf (Basic reduction, image stacking, Image align)  
Alipy (Image align)  

**Other softwares :**  
SExtractor (E.Bertin 1996)  
(Scamp)  
Astrometry.net  
Hotpants  
wcsremap  

Install
-------------  

Contact
-------------
Please don't hesitate to send me an email --> lim9gu@gmail.com
(Suggestions, Questions, Complains, Compliments, Collaboration, Funding...)

Acknowledgement
------------- 
Superviser : Myunshin Im (CEOU, SNU)  
Technical advices : Jinguk Seo (SNU), Kihyun Kim (SNU)  
IMSNG team : Changsu Choi (CEOU, SNU), Gregory S.H. Paek (CEOU, SNU)  
Software development : Yoonsoo Bach Park (SNU)
