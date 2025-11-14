README.txt

PROJECT 2:

link to online github repository: https://github.com/gurmehersk/astr142_proj2


Welcome to the everlong project 2 code! There are a lot of references and notes I want to talk about here!


**UPdate to part a, i changed the fits files and they run better now, no greenness. However, the resolution in pdfs have dropped a bit. I presume the other ones were of a higher quality.

a) The rgb composite might be a little green, displeasing to the eyes: I am aware, and I apologize! However, due to the slow speed of my laptop, which was taking a while to download the 6GB files, and after downloading a wrong file once [i realized after waiting for 30 minutes], I decided to use files which for which shapes could be easily tested and image stacks could be created immediately. The three files i found were found in the link below, and were in the F435W, F606W, and F775W bands. Clearly, this range makes the images a little greener than what would be comfortable. However, it is clear still that some sort of stacking has occured because we can still see the redness in the rgb composite. Admittedly, the blue part of the spectrum is difficult to see, even though I had F435W. I couldn't really change this as other files were giving me errors for some reason, or were taking too long to download. [You will see, i tried multiple different files that have been commented out in the main section of the code]


b) You may notice a lot of code redundancy and code repetition. This occured mainly because I did this project on separate days, and frankly, due to poor documentation initially, I just forgot when I had done some things.

Additionally, it is true that in some cases, even if we were creating a similar figure, we weren't creating the exact same figure! For instance, the overlay plot with the galaxies required different labels, headers, titles than the normal composite image. However, we used the rgb_image function to create the image always, so in that way, there was less code redundancy. 

Could I have probably used other parts of the code and called other functions that returned figures in my code? Yes, i think, maybe? Did I? Not always, just because I was unsure whether it would work. At that point, it was easier to undergo a little bit of redundancy to be 100% sure that the code would work.

c) Throughout the project, I use functions that we mainly saw in class, so there are no surprises. However, whenever I used a function that we had not seen in class, I added the reference link or where I got it from within the comments itself. A few examples of functions i used which weren't introduced: 

- lupton_rgb from astropy (https://docs.astropy.org/en/stable/visualization/rgb.html)
- wcs's pixel_scale_matrix function 
- skycoord_to_pixel [not sure if we saw this in class]
- np.ma.fill --> to clean masked values in our z values 


d) For the cross-verification of catalogs, I discussed various modules like the np.core.records.fromarrays, and the xmatch module from Astropy. However, I thought of an intuitive way that seemed to work best for me. I explained this in detail in the documentation but I used 
.match_to_catalog_sky from skycoords to get the nearest coordinate match from the photometry catalog to the spectroscopy catalog. Before doing this, however, I realized that the number of significant figures in the ra and dec of both catalogs were different, so i rounded the catalogs to 6 decimal places for consistency. Now, if there was a match, the value of the .match_to_catalog_sky would honestly just be 0, since they seemed to perfectly match for a few points that i tested for eye. Still, i gave a 0.0036 arcsecond buffer in the catalog matching, and did the matching using a rather intricate method using np.where(mask) and a basic mask which returned true if the separation < 0.0036. After this was done, I saved the indices where the Mask gave me a true value in the spectrscopy catalog and cross matched those indices with the photometry catalog: I.e., i took advantage of the fact that all the photometric redshift galaxies were present in the spectroscopy catalog, but not vice versa. 


e) For the insets! The most challenging part of this project! I wasn't able to plot the coordinates next to the inset as the wcs scaling changes once we extract the subregion because the pixels now make up different amounts of the RA and DEC. I initially tried creating a new WcS, copying the old one and trimming it somehow, but this just was NOT working. I decided to just plot the pixels as is and remove the tic labels to avoid any irrevelant and redundant pixel axes.  

References:

Fits files: https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
Spectroscopy catalog: https://www.aanda.org/articles/aa/full_html/2017/12/aa31195-17/aa31195-17.html (Inami+2017)
Photometry catalog : https://ui.adsabs.harvard.edu/abs/2016yCat..51500031R/abstract (Rafelski+2015)
np.ma.filled : https://numpy.org/doc/2.3/reference/generated/numpy.ma.filled.html