from PIL import Image

img = Image.open("RamanSpectrumAnalysis_icon.png").convert("RGBA")

img.save(
    "RamanSpectrumAnalysis.ico",
    format="ICO",
    sizes=[
        (16, 16),
        (24, 24),
        (32, 32),
        (48, 48),
        (64, 64),
        (128, 128),
        (256, 256),
    ]
)