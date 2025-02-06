{
  description = "Vlasiator developement environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs, ... }@inputs:
    let 
      pkgs = import nixpkgs { system = "x86_64-linux"; };
      vlsvSrc = builtins.fetchGit {
        url = "https://github.com/fmihpc/vlsv";
        ref = "master";
      };
      
      phiprofSrc = builtins.fetchGit {
        url = "https://github.com/fmihpc/phiprof";
        ref = "master";
      };
      
      zoltanSrc = builtins.fetchGit {
        url= "https://github.com/trilinos/Trilinos.git";
        ref = "master";
      };

      papiSrc = builtins.fetchGit {
        url = "https://github.com/icl-utk-edu/papi.git";
        ref = "master";
      };
     
      vlsvPkg = pkgs.stdenv.mkDerivation {
        pname = "vlsv";
        version = "latest";

        src = vlsvSrc;

        nativeBuildInputs = [ pkgs.makeWrapper ];
        buildInputs = [ pkgs.gcc pkgs.openmpi ];

        buildPhase = ''
          make -j 4
        '';

        installPhase = ''
          mkdir -p $out/lib
          mkdir -p $out/include
          cp -r ./* $out/lib/
          cp -r ./* $out/include/
        '';
      };
      
      phiprofPkg = pkgs.stdenv.mkDerivation {
        pname = "phiprof";
        version = "latest";

        src = phiprofSrc;

        nativeBuildInputs = [ pkgs.makeWrapper ];
        buildInputs = [ pkgs.gcc pkgs.openmpi ];

        buildPhase = ''
          cd src/
          make -j 4
        '';

        installPhase = ''
          cd ..
          mkdir -p $out/lib $out/include
          cp  lib/* $out/lib/
          cp  include/* $out/include/
        '';
      };
      
      zoltanPkg = pkgs.stdenv.mkDerivation {
        pname = "zoltan";
        version = "latest";

        src = zoltanSrc;

        nativeBuildInputs = [ pkgs.cmake];
        buildInputs = [ pkgs.gcc pkgs.openmpi pkgs.perl ];


        configurePhase = ''
          mkdir zoltan_build
          cd zoltan_build
          cmake .. -DCMAKE_INSTALL_PREFIX:FILEPATH=$out -DTrilinos_ENABLE_Zoltan:BOOL=ON    -DZoltan_ENABLE_ULLONG_IDS:Bool=ON  -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF  -DTPL_ENABLE_MPI:BOOL=ON
        '';

        buildPhase = ''
          make -j 8
        '';

        installPhase = ''
          make install
        '';
      };

       # Fetch and build jemalloc
      jemallocPkg = pkgs.stdenv.mkDerivation {
        pname = "jemalloc";
        version = "5.3.0";

        src = pkgs.fetchurl {
          url = "https://github.com/jemalloc/jemalloc/releases/download/5.3.0/jemalloc-5.3.0.tar.bz2";
          sha256 = "sha256-LbgtHnEZ3z5xt2QCGbbf6EeJvAU3mDw7esT3GJrs/qo="; # Update this hash if necessary
        };

        nativeBuildInputs = [ pkgs.gnumake pkgs.autoconf pkgs.automake pkgs.libtool ];
        buildInputs = [ pkgs.gcc ];

        configurePhase = ''
          tar xf $src
          cd jemalloc-5.3.0
          ./configure --prefix=$out --with-jemalloc-prefix=je_
        '';

        buildPhase = ''
          # cd jemalloc-5.3.0
          make -j 8
        '';

        installPhase = ''
          # cd jemalloc-5.3.0
          make install
        '';
      };

      papiPkg = pkgs.stdenv.mkDerivation {
        pname = "papi";
        version = "latest";

        src = papiSrc;

        nativeBuildInputs = [ pkgs.gnumake pkgs.autoconf pkgs.automake pkgs.libtool ];
        buildInputs = [ pkgs.gcc ];

        configurePhase = ''
          cd src
          ./configure --prefix=$out
        '';

        buildPhase = ''
          make -j 8
        '';

        installPhase = ''
          make install
        '';
      };



      boostPkg = pkgs.stdenv.mkDerivation {
        pname = "boost";
        version = "1.72.0";

        src = pkgs.fetchurl {
          url = "http://freefr.dl.sourceforge.net/project/boost/boost/1.72.0/boost_1_72_0.tar.bz2";
          sha256 = "sha256-WcmydLxFHPkam6HdLH/cr11gsbOqg/LJ+hQ0F8xmByI="; # Update this hash if necessary
        };

        nativeBuildInputs = [ pkgs.gnumake pkgs.perl pkgs.bash ];
        buildInputs = [ pkgs.gcc pkgs.openmpi ];

        configurePhase = ''
          tar xf $src
          cd boost_1_72_0
          ./bootstrap.sh --with-libraries=program_options
          echo "using mpi ;" >> ./tools/build/src/user-config.jam
        '';

        buildPhase = ''
          ./b2
        '';

        installPhase = ''
          ./b2 --prefix=$out install
        '';
      };



    in
    {
      devShells.x86_64-linux.default = pkgs.mkShell {
        name = "dev-shell";
        buildInputs = [
          pkgs.git
          pkgs.gcc
          pkgs.openmpi
          vlsvPkg
          phiprofPkg
          zoltanPkg
          jemallocPkg
          papiPkg
          boostPkg
        ];
        shellHook = ''
           export LIB_VLSV=-L${vlsvPkg}/lib
           export INC_VLSV=-I${vlsvPkg}/include
           export LIB_PROFILE=-L${phiprofPkg}/lib/
           export INC_PROFILE=-I${phiprofPkg}/include/
           export LIB_ZOLTAN=-L${zoltanPkg}/lib
           export INC_ZOLTAN=-I${zoltanPkg}/include
           export LIB_JEMALLOC=-L${jemallocPkg}/lib
           export INC_JEMALLOC=-I${jemallocPkg}/include
           export LIB_PAPI=-L${papiPkg}/lib
           export INC_PAPI=-I${papiPkg}/include
           export LIB_BOOST=-L${boostPkg}/lib
           export INC_BOOST=-I${boostPkg}/include
           export VLASIATOR_ARCH=nix
           echo "Vlasiator enviroment ready! Use make -j <cores> to build a fresh vlasiator binary!"
        '';

      };
    };
}


