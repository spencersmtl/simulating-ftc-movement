# Run this to save the output of eigenpatch-gp-landscape. Saves as true number of pixels without caption
ggsave("save.png", patch_map + theme(plot.caption = element_blank()), width = 100/300, height = 100/300, dpi = 300)
# Run this to save the pixel and patch information
saveRDS(gp_landscape, "save_data.rds")
