import qupath.ext.stardist.StarDist2D
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.scripting.QP
import qupath.lib.images.servers.ImageServer
import qupath.lib.images.servers.bioformats.BioFormatsImageServer
import qupath.lib.images.servers.TransformedServerBuilder
import qupath.lib.gui.viewer.QuPathViewer
import java.awt.image.BufferedImage
import java.io.File
import javax.imageio.ImageIO
import qupath.lib.objects.PathAnnotationObject
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

def pathObjects = QP.getSelectedObjects()

def image = getCurrentImageData()
def imageServer = getCurrentServer()
def list = ['DAPI','PD-L1','CD66b','MIF','APOC3','MECA79','CD20','S100A12','Podoplanin','SMA','CD45',
    'CD21','CD31','CD68','CD8','S100A8_9','CD4','MPO','Pan-Cytokeratin','CD3e','ISG15','HLA-DR','CD74']

for (int j = 0;j<100;j++) {
    Annotation = pathObjects[j]
    print Annotation
    def folderPath = "E:\\01-TLS\\09-CODEX Alignment\\04-Qupath Out\\01-135121\\" + Annotation
    Path folder = Paths.get(folderPath)
    Files.createDirectories(folder)
    for (int i = 0;i<23;i++) {
        channel = list[i]
        def outputPath = folderPath+'\\'+channel+ '.tif'
        def roi = Annotation.getROI()
        def singleChannel = new TransformedServerBuilder(image.getServer())
                .extractChannels(channel)
                .build()
        def regionRequest = RegionRequest.createInstance(singleChannel.getPath(),1, roi)
        def roiImage = singleChannel.readRegion(regionRequest)
        def outputFile = new File(outputPath)
        writeImage(roiImage,outputPath)
        println("ROI图像已保存为TIFF文件：$outputPath")
}

}