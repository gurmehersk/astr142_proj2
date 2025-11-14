README.txt

PROJECT 2:

Welcome to the everlong project 2 code! There are a lot of references and notes I want to talk about here!

a) The rgb composite might be a little green, displeasing to the eyes: I am aware, and I apologize! However, due to the slow speed of my laptop, which was taking a while to download the 6GB files, and after downloading a wrong file once [i realized after waiting for 30 minutes], I decided to use files which for which shapes could be easily tested and image stacks could be created immediately. The three files i found were found in the link below, and were in the F435W, F606W, and F775W bands. Clearly, this range makes the images a little greener than what would be comfortable. However, it is clear still that some sort of stacking has occured because we can still see the redness in the rgb composite. Admittedly, the blue part of the spectrum is difficult to see, even though i had F435W. I couldn't really change this as other files were giving me errors for some reason, or were taking too long to download.

b) You may notice a lot of code redundancy and code repetition. This occured mainly because I did this project on separate days, and frankly, due to poor documentation initially, I just forgot when I had done some things.

Additionally, it is true that in some cases, even if we were creating a similar figure, we weren't creating the exact same figure! For instance, the overlay plot with the galaxies required different labels, headers, titles than the normal composite image. However, we used the rgb_image function to create the image always, so in that way, there was less code redundancy. 

Could I have probably used other parts of the code and called other functions that returned figures in my code? Yes, i think, maybe? Did I? Not always, just because I was unsure whether it would work. At that point, it was easier to undergo a little bit of redundancy to be 100% sure that the code would work.

c) Throughout the project, I use functions that we mainly saw in class, so there are no surprises. However, whenever I used a function that we had not seen in class, I added the reference link or where I got it from within the comments itself. A few examples of functions i used which weren't introduced: 

- lupton_rgb from astropy (https://docs.astropy.org/en/stable/visualization/rgb.html)
- wcs's pixel_scale_matrix function 
- skycoord_to_pixel [not sure if we saw this in class]
- np.ma.fill --> to clean masked values in our z values 

References:

Fits files: https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
MUSE Spectrograph catalog: https://www.aanda.org/articles/aa/full_html/2017/12/aa31195-17/aa31195-17.html (Inami+2017)
Photometry catalog : https://ui.adsabs.harvard.edu/abs/2016yCat..51500031R/abstract (Rafelski+2015)
