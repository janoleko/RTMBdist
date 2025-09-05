library(ggplot2)
library(hexSticker)
library(RTMBdist)
library(scales)
library(colorspace)

x <- seq(-1, 4, length = 300)
d <- dskewt(x, 0, 1, 4, 4)
df <- data.frame(x = x, d = d)

p <- ggplot(df, aes(x, d)) +
  geom_area(fill = alpha("deepskyblue", 0.3)) +
  geom_line(linewidth = 1, color = "deepskyblue") +
  theme_void()

sticker(
  p,
  package = "RTMBdist",
  s_x = 1, s_y = 1.21,
  s_width = 1.3, s_height = 1,
  p_size = 20,          # font size
  p_color = alpha("black", 0.8), # text color
  p_y = 0.56,            # vertical position (0 = bottom, 1 = top)
  p_family = "sans",    # font family
  p_fontface = "bold",
  h_fill = lighten("deepskyblue", 0.8),
  h_color = "deepskyblue",
  filename = "./man/figures/RTMBdist_hex.png"
)

