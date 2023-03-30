Games101 课程实验进行中

- HW1, 比较easy
- HW2, z-buffer光栅化，MSAA，MSAA实现得不太好，偷懒改了framework里depth_buf的类型。HW2的代码框架存在一些错误，写的时候还没有意识到，所以代码漏洞可能较多，这部分都会在HW3中修正。
- HW3, 模拟完整的渲染pipeline，记录几个关键点
  - 模型load成triangle_list，然后开始对每个triangle执行draw。
  - draw的过程：
    1. 三角形顶点的mvp变换->齐次除法->视口变换。这部分细节在于顶点在缩放过程可能会导致法线偏移，所以要用专门的一个矩阵去处理法向量，在将法向量赋给视口变换后的三角形。此外对视口变换的理解也很重要，进行齐次除法后，(x,y,z)投影到[-1,1]的NDC空间，而w中保存的就是原3d点的深度值。视口变换把NDC空间的(x,y)映射到图片像素点，z则是z=z*f1+f2，其中f1=(far-near)/2,f2=(near+far)/2,看了一些博客说是做齐次裁剪用，这里用不到。所以视口变换后的点(x,y,z,w)，原深度一直是存在w里的。（我一开始还以为是在z里，拿去插值就错了）。
    2. 我们获得了三角形的顶点在像素坐标系上的位置，接下来就可以光栅化了。使用像素坐标在三角形上的重心坐标来插值三角形顶点的各种属性（坐标、法线、纹理坐标）。你可能会想，我的三角形明明是三维空间的，我只有二维空间的像素坐标，怎么插值呢？首先算重心坐标，这个是在二维空间上算出来的，使用三维三角形的x，y维度。然后使用存在w里的深度计算透视矫正插值，就可以用像素面上投影的二维重心坐标来插值三维三角形的属性了，具体代码去看作业3的光栅化函数。（原games101框架里的透视矫正写错了，这版代码里做了修改）
    3. 作业实现的其实是phong shading，也就是针对像素单位的着色，所以我们是在三维空间里插值了三角形的法向量作为像素的法向量。在shading模型中，三角形只经过mv变换后的（透视投影前的）顶点插值被传入作为shading point，光源（位置，强度）手写定义，视点就是相机位置（框架手写定义的相机位置有小问题）。ok，以上各部分都有了，对着bulin-phong模型的图写就行了。
