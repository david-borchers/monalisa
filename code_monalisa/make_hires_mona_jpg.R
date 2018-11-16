## Making the smaller 1200 x 1200 sub-image of the Mona Lisa, as used
## in later analyses

library(imager)

mona <- load.image("monalisa.jpg")
gray_mona <- grayscale(mona)
small_mona <- resize(gray_mona, -50, -50)

small_mona_face <- imsub(small_mona, y < 1300)
plot(small_mona_face)
dim(small_mona_face)
# user sets desired width and height of image 
desired_width <- 1200
desired_height <- 1200

# how many pixels to trim off each side of image
x_pixels_to_trim <- (width(small_mona_face) - desired_width) 
y_pixels_to_trim <- (height(small_mona_face) - desired_height) 

# horizontal crop removes same/similar number of pixels from each side
if(x_pixels_to_trim %% 2 == 0){ # if number of pixels to trim is even, same from both sides
  small_mona_face <- imsub(small_mona_face, 
                           x > (x_pixels_to_trim / 2) & 
                             x <= ((width(small_mona_face) - x_pixels_to_trim / 2)))
} else { # if number of pixels to trim is odd, take one more from one side
  small_mona_face <- imsub(small_mona_face, 
                           x > floor(x_pixels_to_trim / 2) & 
                             x <= ((width(small_mona_face) - ceiling(x_pixels_to_trim / 2))))
}

# vertical crop removes same/similar number of pixels from each side
if(y_pixels_to_trim %% 2 == 0){ # if number of pixels to trim is even, same from both sides
  small_mona_face <- imsub(small_mona_face, 
                           y > (y_pixels_to_trim / 2) & 
                             y <= ((height(small_mona_face) - y_pixels_to_trim / 2)))
} else { # if number of pixels to trim is odd, take one more from one side
  small_mona_face <- imsub(small_mona_face, 
                           y > floor(y_pixels_to_trim / 2) & 
                             y <= ((height(small_mona_face) - ceiling(y_pixels_to_trim / 2))))
}

dim(small_mona_face)
plot(small_mona_face)

save.image(small_mona_face, file = "output/hires_mona.jpg")
