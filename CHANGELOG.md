# Changelog

## [1.2.0](https://github.com/MPUSP/mpusp-snakemake-wrappers/compare/v1.1.0...v1.2.0) (2025-08-01)


### Features

* added fimo search as third module for MEME suite ([4ebb227](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/4ebb227e9af13f38f3d72bbbef23b1e49094e293))
* added logomaker wrapper ([565bb5f](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/565bb5f28399678838642010ad81004b45fa142d))
* added logoplot wrapper ([3d53752](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/3d5375231c99261b4f233cbd5b74c06846dcd4b3))
* added meme module to generate pwm ([edaa975](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/edaa9752fdddb81dbd3e51ce36b2e62119778bb8))
* added meme module to generate pwm, make motif, and search for motif in target sequence ([e248df9](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/e248df9670ec3d669b25b5b7f270e543b7972121))
* added second MEME suite module to create meme motif ([b5a4d5c](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/b5a4d5cef83c4a6957d750edf07af472803925af))


### Bug Fixes

* added missing pin, formatting ([4cf3595](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/4cf359508db44a827b4bacc8c6f5404bbb6658ea))
* missing dependency ([bf80934](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/bf80934f7c38c3515d5525f6cf65436a5b1feb20))

## [1.1.0](https://github.com/MPUSP/mpusp-snakemake-wrappers/compare/v1.0.0...v1.1.0) (2025-07-28)


### Features

* added new test to get genome locally ([85307c9](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/85307c9e418653faff251e94e5ed6d84c1af7928))
* added test ([4cc8fb3](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/4cc8fb3fe779573af9df6ae59491c082edb8869a))


### Bug Fixes

* added missing fasta indexing ([4cc8fb3](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/4cc8fb3fe779573af9df6ae59491c082edb8869a))
* added missing fasta indexing ([f62b3e8](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/f62b3e8bef10619a02ec6d91396923b26057556f))

## 1.0.0 (2025-06-04)


### Features

* added documentation and CI actions ([65698ac](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/65698acee9728564a7e0f09b4740b4aff62ed27b))
* added dwgsim wrapper ([6eb7935](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/6eb7935fd3ab458074eec147b675107ce20b41c7))
* added get_genome wrapper ([0ae26a7](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/0ae26a7333e16f6ea1c6e48ee4407181538a4715))
* added gffcompare wrapper using official wrapper format ([3035c57](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/3035c57873604f9feadffa251140f6a342d39d62))
* added GH action to test workflows automatically ([a30388e](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/a30388e02b3605751a533aaf4faef2d98dbf9cf5))
* added nanoplot wrapper, input handling ([f16b50c](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/f16b50c2197b1ed0650cdac54170a8ee587ad500))
* added new dwgsim wrapper ([3d0d608](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/3d0d608c55fa07edffb51f071d11e8418353761e))
* added new pycoQC wrapper ([b202111](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/b20211199c1729da0f59970a915b8fb75e3b06df))
* added README ([eb2fe6c](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/eb2fe6c25cb8b28dafd913db559df2a67886ae73))
* added release please workflow ([6ca04cf](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/6ca04cfbdef99ec43f76fe906d1359a719d53c28))
* added remaining files for nanoplot wrapper ([40a8355](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/40a8355ef3f6b64796909f2ca157be10f966b9ad))
* added stringtie assemble wrapper, minor fixes ([f6c8e5f](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/f6c8e5f7e772f9b83f7c020abc19d26568251e7a))
* added StringTie wrapper ([8b879a3](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/8b879a3193e28bb5a4d4788ecfa72b841bfc843d))
* new stringtie merge wrapper ([c1e7b1b](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/c1e7b1bd15cda85300bea5a864f8ecdd5b4c5472))
* new wrappers, new fully automatic guithub action test suite ([5fcbd05](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/5fcbd050ce041bb678ad768dadc1a2912df4d4ab))
* update GH actions workflow to create list of wrapper tests automatically ([32e4e41](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/32e4e417bca1328976d0f6bd826fe117de765fb4))


### Bug Fixes

* add statement to renumber features, closes [#7](https://github.com/MPUSP/mpusp-snakemake-wrappers/issues/7) ([1eba57d](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/1eba57d91b4ed3a9468b8832bafb151b572fceaf))
* added 1 optional input and improved input handling ([f3b36a3](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/f3b36a3da7b98b4bfa24629f612e2af64df05700))
* added choice between classic and strand specific ([1f24277](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/1f242775e85bcf9bd144e947b81e6848d5900eeb))
* added missing log file ([33f519b](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/33f519bf1dd10a5e1feb9282f2a005989a3d22c9))
* added missing pinned env ([6e4e2d8](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/6e4e2d803f8bf2c7be6c93d4028c50a2e750defa))
* added more info on test and test result badge ([b036fc0](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/b036fc00b57ae5cd693647971856e3344c1f9117))
* added threads ([1cea51a](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/1cea51a4d9452274d0451cc3a017dacae8c6a73c))
* another try to enable snakemake CI test ([b7983bc](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/b7983bc9c672565132c460879bf2feaacb6ca5f3))
* documentation to suppress additional report files ([faa1992](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/faa19928f8cffce197ecd9b7a9dc8f91c667d371))
* fixed sorting to have tx first ([049851e](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/049851e9bf3d0b48aa79b067301e1c8ce49a8d34))
* formatting ([19bad4a](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/19bad4ae576a1560248f11391703944d0c405d76))
* formatting single input files ([ae14271](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/ae14271bfa43d04f9acc255ca81d9d3d08340bbc))
* input and log path ([6bce915](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/6bce91587b7886be241389f9b47f3b4271957d5a))
* input and log path v2 ([cd79fa9](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/cd79fa9835a19fde591827fabfc516e85f25505d))
* linting ([e43e97d](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/e43e97dc071face1d595c2e7ee88f8c954f45ce6))
* linting error ([63be326](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/63be32630926145d5929b146d24356c3f6fb1916))
* made stringtie merge strand specific ([c456ecf](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/c456ecfdfaa7c5bfeea9a794fcd5b1728f503fb2))
* made stringtie merge strand specific ([04b2601](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/04b2601d1a2210f3c4e3230ff089d7ba3a95faad))
* minor improvements to usability ([d076c1f](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/d076c1f69df508b3e18b47128f82e2d1af44f9ff))
* missing param docs ([de8b1c2](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/de8b1c2fa17fa4c985fde94142f0f9231291477f))
* more debugging of actions workflow ([4c879ad](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/4c879adfdc0d0c21005aad0ca4e01e68021eb936))
* release workflow, bug fix ([8987068](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/89870681b441e5a615ac0f488d9ce58d3580da6a))
* renumbering and sorting stringtie merge gff output, and new pycoQC wrapper ([4c3eb9e](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/4c3eb9edf5bb3d06ea33e3da620edf85116aa525))
* unzipping error based on wrong cli args ([b053999](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/b0539999dcd1a1e1ec51c9c3000e3d2c949266c5))
* unzipping error based on wrong cli args ([60b0dc8](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/60b0dc8c000b1a1e371823cc90cf8015c121a213))
* update GH actions artifact ([21229b0](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/21229b02cf32e8ee652a79f01b83f111dbe0323a))
* updated GH action tests ([883b777](https://github.com/MPUSP/mpusp-snakemake-wrappers/commit/883b777cbfbd85e61b6c0f646dae01ecf233e660))
