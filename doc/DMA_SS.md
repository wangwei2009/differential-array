## 双麦差分+谱减降噪

[TOC]

### 1.双麦配置[^1]

<img src="E:\work\matlab\Github\elite\2MIC_NS\pcm\differential-array\doc\pic\57563270-4C95-4f4e-A43E-A811F6E1ECE2.png" alt="57563270-4C95-4f4e-A43E-A811F6E1ECE2" style="zoom:50%;" />

### 2.处理流程[^2]

<img src="E:\work\matlab\Github\elite\2MIC_NS\pcm\differential-array\doc\pic\{97F38125-6E72-4774-B031-BF57999C8F17}_20191021153136.jpg" alt="{97F38125-6E72-4774-B031-BF57999C8F17}_20191021153136" style="zoom:50%;" />

处理流程分为固定波束和谱减两个部分

#### 2.1 固定波束

已假定目标方向为正前方

$b_{12}$、$b_{21}$和$n_{12}$分别表示消掉延时为$+\tau$、$-\tau$的信号，表示消掉目标信号后的参考噪声

​    1. $n_{12}$  时域公式表示为：
$$
n_{12}(t)=x_{1}(t)-x_{2}(t)\tag1
$$


​       频域表达式为：
$$
Y_{n_{12}}(\omega)=\left[\begin{array}{c}{1} \\ {-1}\end{array}\right] S(\omega)\tag2
$$




2. $b_{12}$和$b_{21}$分别用来消掉从端射方向入射的信号

   时域表达式为：

$$
b_{12}(t)=x_{1}(t-\tau)-x_{2}(t)\\
b_{21}(t)=x_{1}(t)-x_{2}(t-\tau)\tag3
$$

​       频域表达式为：
$$
Y_{b_{12}}(\omega)=\left[\begin{array}{c}{e^{-j \omega \tau}} \\ {-1}\end{array}\right] S(\omega)\\
Y_{b_{21}}(\omega, \theta)=\left[\begin{array}{c}{e^{-j \omega \tau}} \\ {-1}\end{array}\right] S(\omega)\tag4
$$

3. 波束输出

信号从$\theta$角入射，则上面三个固定波束输出分别表示为
$$
N_{12}(\omega, \theta)=\left[\begin{array}{ll}{1} & {e^{-j \omega \tau_{0} \sin \theta}}\end{array}\right]\left[\begin{array}{c}{1} \\ {-1}\end{array}\right] S(\omega)\tag5
$$

$$
B_{12}(\omega, \theta)=\left[\begin{array}{ll}{1} & {e^{-j \omega \tau_{0} \sin \theta}}\end{array}\right]\left[\begin{array}{c}{e^{-j \omega \tau}} \\ {-1}\end{array}\right] S(\omega)\tag6
$$

$$
B_{12}(\omega, \theta)=\left[\begin{array}{ll}{1} & {e^{-j \omega \tau_{0} \sin \theta}}\end{array}\right]\left[\begin{array}{c}{1} \\ {-e^{-j \omega \tau}}\end{array}\right] S(\omega)\tag7
$$

画出三个固定波束图如下

<img src="E:\work\matlab\Github\elite\2MIC_NS\pcm\differential-array\doc\pic\{73B28CF7-E00C-46C0-8393-FFED146AAD1D}_20191021203234.jpg" alt="{73B28CF7-E00C-46C0-8393-FFED146AAD1D}_20191021203234" style="zoom:50%;" />

#### 2.2 谱减

用$M_{12}$对目标方向加强
$$
\left|M_{12}(\omega, k)\right|=\min \left[\left|B_{12}(\omega, k)\right|, | B_{21}(\omega, k)\right|]\tag8
$$
$M_{12}$的输出如下所示

<img src="E:\work\matlab\Github\elite\2MIC_NS\pcm\differential-array\doc\pic\{EC4ED9E9-1EC2-4136-BB57-76AF5EE9E4F3}_20191023141233.jpg" alt="{EC4ED9E9-1EC2-4136-BB57-76AF5EE9E4F3}_20191023141233" style="zoom:50%;" />

功率谱减公式如下：
$$
\left|Y_{12}^{\prime}(\omega, k)\right|^{2}=\left\{\begin{array}{ll}{\left|M_{12}(\omega, k)\right|^{2}-\left|N_{12}(\omega, k)\right|^{2},} & {\text { if }\left|M_{12}(\omega, k)\right|>\left|N_{12}(\omega, k)\right|} \\ {0,} & {\text { otherwise }}\end{array}\right.\tag9
$$
到这一步的输出波束图如下：

<img src="E:\work\matlab\Github\elite\2MIC_NS\pcm\differential-array\doc\pic\{8CEB0A6C-7635-45BE-9668-B327E884334C}_20191023141623.jpg" alt="{8CEB0A6C-7635-45BE-9668-B327E884334C}_20191023141623" style="zoom:50%;" />

#### 2.3 频率补偿

从$M_{12}$的波束图可以看到，当信号从正前方入射时，输出响应并不为1，最后需要添加一个补偿滤波器，分析过程参考[1]，补偿系数形式如下
$$
H_L(\omega)=\frac{1}{\sqrt{2(1-\cos{\omega\tau_0})}}\tag{10}
$$
经过频率补偿后，最终的频率不变波束图如下

<img src="E:\work\matlab\Github\elite\2MIC_NS\pcm\differential-array\doc\pic\{5FF73B88-7BB7-4AE4-9EC0-37BBDD2499F8}_20191023145635.jpg" alt="{5FF73B88-7BB7-4AE4-9EC0-37BBDD2499F8}_20191023145635" style="zoom:50%;" />

#### 2.4  相位恢复

和常规谱减法类似，使用带噪信号的相位作为目标信号的相位

### 3.误差分析

以上步骤最终画出来的波束图很理想，但实际录音数据测试效果并不好，分析可能原因如下

* DOA误差
* 混响影响
* 双麦间距
* 谱减

#### 3.1 DOA误差

从最终的波束图可以看到，目标方向非常窄，用仿真数据测试，在$\pm20$度的时候目标声音就有一定程度的衰减变小

#### 3.2 混响影响

用仿真数据测试，当目标声音没有混响时，输出质量很好，几乎听不到失真，当混响阶数及混响时间加大时，输出失真相应增大，分析可能原因为当有混响时，输出的参考噪声$N_{12}$的还是能听到目标声音，在下一步的谱减中会消掉目标声音而造成失真。

#### 3.3 双麦间距

差分阵列对间距有要求，间距过大时无法画出频率不变波束图，仿真测试0.005间距的效果比0.025间距的失真小

#### 3.4 谱减法

(9)式是基础的半波整形谱减法，容易产生音乐残留噪声，可以尝试使用参数控制的过减谱减以及多带谱减等改进

### 4. 改进

贴一张有干扰录音处理前后的频谱图

<img src="E:\work\matlab\Github\elite\2MIC_NS\pcm\differential-array\doc\pic\{1D9ED7CA-80A0-42BE-A554-B9F6A638282F}_20191023164007.jpg" alt="{1D9ED7CA-80A0-42BE-A554-B9F6A638282F}_20191023164007" style="zoom:33%;" />



<img src="E:\work\matlab\Github\elite\2MIC_NS\pcm\differential-array\doc\pic\{E5D9A13E-3F78-42D2-B864-DF4168E58AC6}_20191023164103.jpg" alt="{E5D9A13E-3F78-42D2-B864-DF4168E58AC6}_20191023164103" style="zoom:33%;" />

​		批处理实际录音，识别效果并没有明显提升，根据以上分析，可以尝试的改进方法为在前一级加去混响模块、改进谱减算法减小音乐噪声等。

​		目前仿真数据的效果还是要好于实际录音的效果，可能还有其它因素没有考虑到，算法计算过程可能还需要深入理解，后续有发现再改进







参考：

[^1]:《Sound Source Separation Using Null-Beamforming and Spectral Subtraction for mobile Device》
[^2]:《differential microphone arrays》