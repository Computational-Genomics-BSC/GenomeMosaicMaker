# CHANGELOG



## v0.1.0 (2023-11-28)

### Feature

* feat: Add split-read-groups argument ([`412c0e6`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/412c0e664b6fcd2147b8c509b75c889a2e799141))

* feat: Preserve RG from reads ([`701d328`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/701d328a052ee493881993c0d4a9ed3d787ac247))

* feat: Add infill file ([`8db3f75`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/8db3f75f96d2913f4a5ba24c42ef7b8dff59fad4))

* feat: Complete multitasking ([`8b9a9bb`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/8b9a9bb7e56a791ec1b8d80b7db2ce42efd94568))

* feat: Add multithreading capabilities ([`1f8b11a`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/1f8b11a76ca26bed86a4346193989c24d388ee68))

* feat: Add source code ([`0044210`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/004421094844d00e9fc6309ba54b691c57372e0a))

* feat: Add dependencies and .gitignore ([`a25eab8`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/a25eab810f1ded8e711f89a1e447105e7b0b588f))

### Fix

* fix: Fix Dockerfile ([`49c48a3`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/49c48a3a24cd7ac4fa3a9fa3a67fc3c42ca03d01))

* fix: Fix infill_read typo ([`fd24ec2`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/fd24ec2177761fe309d82561496df8f48e3182a0))

* fix: Add compatibility with different canvas contigs ([`fb07222`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/fb07222f6b0f89c663bd1830fadb4fb5ba2b65b8))

* fix: Set write mode based on file extension ([`0f18dbc`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/0f18dbcd57de92afe9d33683518bef741687e038))

* fix: Fix str typo ([`db723e6`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/db723e6ba2d5abc0daf74be7e3e9a71c4792ad24))

* fix: Add all reads in the zones ([`880059d`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/880059dddc54d13484c84c080bc98910edae78b0))

* fix: Remove temp files at the end of processing ([`ff29437`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/ff29437da27a1ff8d049c1bda1a73d40550c0f67))

* fix: Get all suplementary/secondary reads ([`d8ec13d`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/d8ec13d9fa79f14eb3ce9fbfcd94edd9c0a29361))

* fix: Select primary, secondary and supplementary ([`ef910b8`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/ef910b87359ba5f2a4e6f9243c6239a215bd75de))

* fix: Fix variant input type ([`b7bf6c6`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/b7bf6c61298a99f806b000bb08d511c2d39a6285))

* fix: Add suplementary reads ([`691ddee`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/691ddee283337c04d2bbb241a8ee20c7e6bb54d2))

* fix: Fix missing headers from infill file ([`5a93a9d`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/5a93a9df332ea3938385fbe23e36da76c2303faa))

* fix: Add missing read groups ([`9fca220`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/9fca220b7b6abb014eaf6a0b139c8c8ecbf3d2dc))

* fix: Combine read groups ([`665ba28`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/665ba28f2c7c8c4779f4ab4bf43bb9628376d0b1))

* fix: Avoid checking temporal SAM files ([`7e8948b`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/7e8948ba30c1ce58f3ac3b041a88206e6f2206dd))

* fix: Ignore reads missing mates ([`b7976b1`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/b7976b1878825453d05f4c1f4932c4922869271b))

* fix: Add RG tags ([`76e1a63`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/76e1a63be541b057ac984d74f4556d7ae8daa3b7))

* fix: Fix type paddings ([`9d72436`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/9d724366e899f1eb0f31f319f36f1cc8b81bd7be))

### Performance

* perf: Improve output I/O performance ([`b244da0`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/b244da08f27f57c90e4dbb77d8012a5eb80e7060))

* perf: Improve mate finding performance ([`5faa43e`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/5faa43e231633cf0aef39f932fcc859844dda2dc))

* perf: Improve performance by splitting the job in different processes ([`645b11b`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/645b11b4927b2ea55c994f2cfbf5316d77d77880))

### Unknown

* deps: Remove variant-extractor dependency ([`b9f2dcb`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/b9f2dcb620eac83c20130cd1bed275c30d9f0a2b))

* Merge branch &#39;sorted&#39; into &#39;main&#39;

perf: Improve performance by splitting the job in different processes

See merge request rmarti1/genome-combinator!1 ([`5aebbaf`](https://github.com/Computational-Genomics-BSC/GenomeMosaicMaker/commit/5aebbaf3aa3aa59a3991f7f7db44a0c6b06b008e))
