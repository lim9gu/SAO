SAO1-m
=============
Data reduction package for Seoul National University Astronomical Observatory (SAO) 1-m telescope. Firstly, this code contains basic reduction of bias/dark subtraction, and flat fielding (skyflat, domeflat). Data classification (bias, Dark, domeflat, skyflat) can be performed. (Ex. Exposure time, Filter, Binning). Solving WCS is based on Astrometry.net software (Using 'Scamp' will be added soon). python2.7 and python3.6 is compatible.   

Current version is 0.5. This code will be updated continuously.   
(Last update: 2019-05-23)  

Prerequisites
-------------  
**Python libraries :**  
astropy  
numpy  
Pyraf (Basic reduction, image stacking, Image align)  
Alipy (Image align)  

**Other softwares :**  
SExtractor (Photometry, Masking, Sky subtraction ; E.Bertin 1996)  
(Scamp, Image stacking, resampling)    
Astrometry.net (WCS solving)  
Stiff (Color image)  
Hotpants (Image subtraction)  
wcsremap (Image alignment)  

Install
-------------
1. Reduction code : Preprocessing functions are contained. Download SAO_process.py on your site-packages/SAO location with __init__.py file.
2. Download code : Using IP and ID of the data storage server, download data on specific date. This code can be erased without prediction due to security. Download SAO_download.py.  
3. Command code : Using reduction code 1., all process of reduction is performed automatically. Download SAO_command.py 
Continuous updates will be applied.

Contact
-------------
Please don't hesitate to send me an email (lim9gu@gmail.com)  
(Suggestions, Questions, Complains, Compliments, Collaboration, Funding...)

Acknowledgement
------------- 
Superviser : Professor. Myungshin Im (CEOU, SNU)  
Technical advices : Jinguk Seo (SNU), Kihyun Kim (SNU)  
IMSNG team : Changsu Choi (CEOU, SNU), Gregory S.H. Paek (CEOU, SNU), Sophia Kim (CEOU, SNU), Joonho Kim (CEOU, SNU)  
Software development : Yoonsoo P. Bach (SNU)
