
import win32com.client

def make_rgb(r, g, b):
    color = win32com.client.Dispatch("Illustrator.RGBColor")
    color.Red = r
    color.Green = g
    color.Blue = b
    return color

ai = win32com.client.Dispatch("Illustrator.Application")

# 新建 A4 文档 (RGB, 宽595pt, 高842pt)
doc = ai.Documents.Add(2, 595, 842)

layer = doc.Layers[0]
layer.Name = "main"

# 画矩形 (top, left, width, height)
rect = layer.PathItems.Rectangle(700, 100, 200, 150)
rect.FillColor = make_rgb(255, 80, 80)
rect.StrokeColor = make_rgb(0, 0, 0)
rect.StrokeWidth = 2

# 画椭圆
ellipse = layer.PathItems.Ellipse(450, 100, 150, 150)
ellipse.FillColor = make_rgb(80, 150, 255)
ellipse.StrokeWidth = 0
ellipse.Stroked = False

# 添加文字
tf = doc.TextFrames.Add()
tf.Position = [100, 800]
tf.Contents = "Hello Illustrator"
tf.CharacterAttributes.Size = 36
tf.CharacterAttributes.FillColor = make_rgb(30, 30, 30)

print("图层数:", doc.Layers.Count)
print("图形数:", layer.PathItems.Count)
print("文字数:", doc.TextFrames.Count)
