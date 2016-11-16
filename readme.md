--
Software requirement
-----
gcc 4.4.7 or 4.8.5

CUDA 5.5, 6.0, 6.5, 7.0, 7.5 (tested and works)

--
Hardware
------
GPU: C2050, C2070, K20, K40 (tested)

--
Compile
-----

make

--
Execute
------
Type: "./gpu_ibfs.bin" it will show you what is needed.

You could use the code from "tuple_text_to_binary_csr" folder to convert a edge list (text formated) graph into CSR files (binary), e.g., if you have a text file called "test.txt", it will convert it into "test.txt_beg_pos.bin" and "test.txt_csr.bin". You will need these two files to run gpu_ibfs.

--
Converter: edge tuples to CSR
----
- Compile: make
- To execute: type "./text_to_bin", it will show you what is needed
- Basically, you could download a tuple list file from [snap](https://snap.stanford.edu/data/). Afterwards, you could use this converter to convert the edge list into CSR format. 

**For example**:

- Download https://snap.stanford.edu/data/com-Orkut.html file. **unzip** it. 
- **./text_to_bin.bin soc-orkut.mtx 1 2(the number may change due to how many lines are not edges)**
- You will get *soc-orkut.mtx_beg_pos.bin* and *soc-orkut.mtx_csr.bin*. 
- You could use these two files to run enterprise.

--
Applications of iBFS
---------
- Centrality computation, e.g., Betweenness centrality, Closeness centrality.
- Multi-source shortest path, All-pairs shortest path.
- Reachiability index construction.
- **In short, any applications, which need to run multiple traversals on the same graph, should benefit**

--
Reference
-------
--
Reference
-------
[SC '15] Enterprise: Breadth-First Graph Traversal on GPUs [[PDF](http://hang-liu.com/publication/enterprise_sc15.pdf)] [[Slides](http://hang-liu.com/publication/enterprise_slides.pdf)] [[Blog](http://hang-liu.com/enterprise_blog.html)]

[SIGMOD '16] iBFS: Concurrent Breadth-First Search on GPUs [[PDF](http://hang-liu.com/publication/ibfs.pdf)] [[Slides](http://hang-liu.com/publication/ibfs_slides.pdf)] [[Poster](http://hang-liu.com/publication/ibfs_poster.pdf)]
