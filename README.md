## Project Information

This project has two main parts.

#### Smith Waterman Algorithm

-  It is a modified Smith Waterman algorithm that can do sequence alignment automatically in python.

#### Introduction

- [Traditional Smith-Waterman Algorithm]([https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm](https://en.wikipedia.org/wiki/Smith–Waterman_algorithm))

It is a dynamic planning algorithms for sequence alignment.

<div align=center>
<img width="190" src="./pictures/image-20200121101831475.png" alt="image-20200121101901465" />

<img width="190" src="./pictures/image-20200121101901466.png" alt="image-20200121101901466" style="zoom:40%;" />
</div>


- Modified Smith-Waterman Algorithm

 We made every two bases in the amino acid as a pair and implemented the new Scoring Matrix of 16*16

##### Output
<div align=center>
<img width="550" src="./pictures/image-20200121101151800.png">
</div>


#### Get Sequence

- It is a crawler in python and bash that can get sequence quickly on NCBI.

##### Output
<div align=center>
<img width="500" src="./pictures/image-20200121102444779.png" alt="image-20200121102444779" style="zoom:50%;" />
</div>
