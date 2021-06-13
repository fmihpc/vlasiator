; ModuleID = 'backgroundfield/quadr.cpp'
target datalayout = "e-p:64:64-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-pc-linux-gnu"
define internal void @pgCplus_compiled.() noinline {
L.entry:
	ret void
}

%struct.T1DFunction = type <{ i32 (...)* (...)*}> 



define internal void @_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a93656trapezERK11T1DFunctionddRdRii(%struct.T1DFunction* %func.arg, double %a.arg, double %b.arg, double* %S.arg, i32* %it.arg, i32 signext %n.arg) #0 !dbg !1359 {
L.entry:
	%func.addr = alloca %struct.T1DFunction*, align 8
	%a.addr = alloca double, align 8
	%b.addr = alloca double, align 8
	%S.addr = alloca double*, align 8
	%it.addr = alloca i32*, align 8
	%n.addr = alloca i32, align 4
	%.Q0000.addr = alloca double, align 8
	%.Q0001.addr = alloca double, align 8
	%delta.addr = alloca double, align 8
	%x.addr = alloca double, align 8
	%sum.addr = alloca double, align 8
	%j.addr = alloca i32, align 4
	%.Q0002.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %func.addr, metadata !1365, metadata !1366), !dbg !1360
	store %struct.T1DFunction* %func.arg, %struct.T1DFunction** %func.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %func.addr, metadata !1367, metadata !1366), !dbg !1360
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1368, metadata !1366), !dbg !1360
	store double %a.arg, double* %a.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1369, metadata !1366), !dbg !1360
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1370, metadata !1366), !dbg !1360
	store double %b.arg, double* %b.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1371, metadata !1366), !dbg !1360
	call void @llvm.dbg.declare (metadata double** %S.addr, metadata !1372, metadata !1366), !dbg !1360
	store double* %S.arg, double** %S.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %S.addr, metadata !1373, metadata !1366), !dbg !1360
	call void @llvm.dbg.declare (metadata i32** %it.addr, metadata !1374, metadata !1366), !dbg !1360
	store i32* %it.arg, i32** %it.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata i32** %it.addr, metadata !1375, metadata !1366), !dbg !1360
	call void @llvm.dbg.declare (metadata i32* %n.addr, metadata !1376, metadata !1366), !dbg !1360
	store i32 %n.arg, i32* %n.addr, align 4, !tbaa !1399
	call void @llvm.dbg.declare (metadata i32* %n.addr, metadata !1377, metadata !1366), !dbg !1360
	%0 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1378
	%1 = icmp ne i32  %0, 1, !dbg !1378
	br i1  %1, label %L.B0000, label %L.B0153, !dbg !1378
L.B0153:
	%2 = load %struct.T1DFunction*, %struct.T1DFunction** %func.addr, align 8, !tbaa !1398, !dbg !1379
	%3 = load double, double* %b.addr, align 8, !tbaa !1403, !dbg !1379
	%4 = bitcast %struct.T1DFunction*  %2 to double (%struct.T1DFunction*, double)***, !dbg !1379
	%5 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %4, align 8, !tbaa !1398, !dbg !1379
	%6 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %5, align 8, !tbaa !1398, !dbg !1379
	%7 = call double  %6 (%struct.T1DFunction*  %2, double  %3), !dbg !1379
	store double  %7, double* %.Q0000.addr, align 8, !tbaa !1403, !dbg !1379
	%8 = load %struct.T1DFunction*, %struct.T1DFunction** %func.addr, align 8, !tbaa !1398, !dbg !1379
	%9 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !1379
	%10 = bitcast %struct.T1DFunction*  %8 to double (%struct.T1DFunction*, double)***, !dbg !1379
	%11 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %10, align 8, !tbaa !1398, !dbg !1379
	%12 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %11, align 8, !tbaa !1398, !dbg !1379
	%13 = call double  %12 (%struct.T1DFunction*  %8, double  %9), !dbg !1379
	store double  %13, double* %.Q0001.addr, align 8, !tbaa !1403, !dbg !1379
	%14 = load double, double* %.Q0000.addr, align 8, !tbaa !1403, !dbg !1379
	%15 = fadd double  %14,  %13, !dbg !1379
	%16 = load double, double* %b.addr, align 8, !tbaa !1403, !dbg !1379
	%17 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !1379
	%18 = fsub double  %16,  %17, !dbg !1379
	%19 = fmul double  %18,  5.00000000000000000E-1, !dbg !1379
	%20 = fmul double  %15,  %19, !dbg !1379
	%21 = load double*, double** %S.addr, align 8, !tbaa !1398, !dbg !1379
	store double  %20, double*  %21, align 8, !tbaa !1399, !dbg !1379
	%22 = load i32*, i32** %it.addr, align 8, !tbaa !1398, !dbg !1380
	store i32 1, i32*  %22, align 4, !tbaa !1399, !dbg !1380
	br label %L.B0001, !dbg !1381
L.B0000:
	%23 = load double, double* %b.addr, align 8, !tbaa !1403, !dbg !1382
	%24 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !1382
	%25 = fsub double  %23,  %24, !dbg !1382
	%26 = load i32*, i32** %it.addr, align 8, !tbaa !1398, !dbg !1382
	%27 = load i32, i32*  %26, align 4, !tbaa !1399, !dbg !1382
	%28 = sitofp i32  %27 to double, !dbg !1382
	%29 = fdiv double  %25,  %28, !dbg !1382
	call void @llvm.dbg.declare (metadata double* %delta.addr, metadata !1383, metadata !1366), !dbg !1364
	store double  %29, double* %delta.addr, align 8, !tbaa !1403, !dbg !1382
	%30 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !1384
	%31 = call double @llvm.fma.f64 (double  %29, double  5.00000000000000000E-1, double  %30), !dbg !1384
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !1385, metadata !1366), !dbg !1364
	store double  %31, double* %x.addr, align 8, !tbaa !1403, !dbg !1384
	call void @llvm.dbg.declare (metadata double* %sum.addr, metadata !1387, metadata !1366), !dbg !1364
	store double  0.00000000000000000E+0, double* %sum.addr, align 8, !tbaa !1403, !dbg !1386
	call void @llvm.dbg.declare (metadata i32* %j.addr, metadata !1389, metadata !1366), !dbg !1360
	store i32 0, i32* %j.addr, align 4, !tbaa !1401, !dbg !1388
	%32 = load i32*, i32** %it.addr, align 8, !tbaa !1398, !dbg !1388
	%33 = load i32, i32*  %32, align 4, !tbaa !1399, !dbg !1388
	%34 = icmp sle i32  %33, 0, !dbg !1388
	br i1  %34, label %L.B0004, label %L.B0003, !dbg !1388
L.B0003:
	%35 = load %struct.T1DFunction*, %struct.T1DFunction** %func.addr, align 8, !tbaa !1398, !dbg !1390
	%36 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1390
	%37 = bitcast %struct.T1DFunction*  %35 to double (%struct.T1DFunction*, double)***, !dbg !1390
	%38 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %37, align 8, !tbaa !1398, !dbg !1390
	%39 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %38, align 8, !tbaa !1398, !dbg !1390
	%40 = call double  %39 (%struct.T1DFunction*  %35, double  %36), !dbg !1390
	store double  %40, double* %.Q0002.addr, align 8, !tbaa !1403, !dbg !1390
	%41 = load double, double* %sum.addr, align 8, !tbaa !1403, !dbg !1390
	%42 = fadd double  %40,  %41, !dbg !1390
	store double  %42, double* %sum.addr, align 8, !tbaa !1403, !dbg !1390
	%43 = load double, double* %delta.addr, align 8, !tbaa !1403, !dbg !1391
	%44 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1391
	%45 = fadd double  %43,  %44, !dbg !1391
	store double  %45, double* %x.addr, align 8, !tbaa !1403, !dbg !1391
	%46 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1392

	%47 = add i32  %46, 1, !dbg !1392
	store i32  %47, i32* %j.addr, align 4, !tbaa !1401, !dbg !1392
	%48 = load i32*, i32** %it.addr, align 8, !tbaa !1398, !dbg !1392
	%49 = load i32, i32*  %48, align 4, !tbaa !1399, !dbg !1392
	%50 = icmp slt i32  %47,  %49, !dbg !1392
	br i1  %50, label %L.B0003, label %L.B0004, !dbg !1392
L.B0004:
	%51 = load double, double* %b.addr, align 8, !tbaa !1403, !dbg !1393
	%52 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !1393
	%53 = fsub double  %51,  %52, !dbg !1393
	%54 = load double, double* %sum.addr, align 8, !tbaa !1403, !dbg !1393
	%55 = fmul double  %53,  %54, !dbg !1393
	%56 = load i32*, i32** %it.addr, align 8, !tbaa !1398, !dbg !1393
	%57 = load i32, i32*  %56, align 4, !tbaa !1399, !dbg !1393
	%58 = sitofp i32  %57 to double, !dbg !1393
	%59 = fdiv double  %55,  %58, !dbg !1393
	%60 = load double*, double** %S.addr, align 8, !tbaa !1398, !dbg !1393
	%61 = load double, double*  %60, align 8, !tbaa !1399, !dbg !1393
	%62 = fadd double  %59,  %61, !dbg !1393
	%63 = fmul double  %62,  5.00000000000000000E-1, !dbg !1393
	store double  %63, double*  %60, align 8, !tbaa !1399, !dbg !1393
	%64 = mul i32  %57, 2, !dbg !1394
	store i32  %64, i32*  %56, align 4, !tbaa !1399, !dbg !1394
	br label %L.B0001
L.B0001:
	ret void, !dbg !1395
}
define internal void @_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a93656polintEPKdS1_idRdS2_(double* %xa.arg, double* %ya.arg, i32 signext %n.arg, double %x.arg, double* %y.arg, double* %dy.arg) #0 !dbg !1408 {
L.entry:
	%xa.addr = alloca double*, align 8
	%ya.addr = alloca double*, align 8
	%n.addr = alloca i32, align 4
	%x.addr = alloca double, align 8
	%y.addr = alloca double*, align 8
	%dy.addr = alloca double*, align 8
	%ns.addr = alloca i32, align 4
	%dif.addr = alloca double, align 8
	%.ndi0002.addr = alloca i32, align 4
	%.lcr011025.addr = alloca i32, align 4
	%dift.addr = alloca double, align 8
	%D.addr = alloca [20 x double], align 8
	%C.addr = alloca [20 x double], align 8
	%.TRP0000.addr = alloca i64, align 8
	%.ndk0004.addr = alloca i64, align 8
	%.lcr051025.addr = alloca i64, align 8
	%.TRP0001.addr = alloca i64, align 8
	%.ndk0006.addr = alloca i64, align 8
	%.lcr051026.addr = alloca i64, align 8
	%.r2.0254.addr = alloca i32, align 4
	%.G0001.addr = alloca i8*, align 8
	%h0.addr = alloca double, align 8
	%hp.addr = alloca double, align 8
	%W.addr = alloca double, align 8
	%den.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata double** %xa.addr, metadata !1418, metadata !1366), !dbg !1409
	store double* %xa.arg, double** %xa.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %xa.addr, metadata !1419, metadata !1366), !dbg !1409
	call void @llvm.dbg.declare (metadata double** %ya.addr, metadata !1420, metadata !1366), !dbg !1409
	store double* %ya.arg, double** %ya.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %ya.addr, metadata !1421, metadata !1366), !dbg !1409
	call void @llvm.dbg.declare (metadata i32* %n.addr, metadata !1422, metadata !1366), !dbg !1409
	store i32 %n.arg, i32* %n.addr, align 4, !tbaa !1399
	call void @llvm.dbg.declare (metadata i32* %n.addr, metadata !1423, metadata !1366), !dbg !1409
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !1424, metadata !1366), !dbg !1409
	store double %x.arg, double* %x.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !1425, metadata !1366), !dbg !1409
	call void @llvm.dbg.declare (metadata double** %y.addr, metadata !1426, metadata !1366), !dbg !1409
	store double* %y.arg, double** %y.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %y.addr, metadata !1427, metadata !1366), !dbg !1409
	call void @llvm.dbg.declare (metadata double** %dy.addr, metadata !1428, metadata !1366), !dbg !1409
	store double* %dy.arg, double** %dy.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %dy.addr, metadata !1429, metadata !1366), !dbg !1409
	call void @llvm.dbg.declare (metadata i32* %ns.addr, metadata !1431, metadata !1366), !dbg !1409
	store i32 0, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1430
	%0 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1432
	%1 = load double*, double** %xa.addr, align 8, !tbaa !1398, !dbg !1432
	%2 = load double, double*  %1, align 8, !tbaa !1399, !dbg !1432
	%3 = fsub double  %0,  %2, !dbg !1432
	%4 = call double @llvm.fabs.f64 (double  %3), !dbg !1432
	call void @llvm.dbg.declare (metadata double* %dif.addr, metadata !1433, metadata !1366), !dbg !1409
	store double  %4, double* %dif.addr, align 8, !tbaa !1403, !dbg !1432
	store i32 0, i32* %.ndi0002.addr, align 4, !tbaa !1401, !dbg !1434
	%5 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1434
	%6 = icmp sle i32  %5, 0, !dbg !1434
	br i1  %6, label %L.B0164, label %L.B0169, !dbg !1434
L.B0169:
	store i32  %5, i32* %.lcr011025.addr, align 4, !tbaa !1401, !dbg !1434
	br label %L.B0163
L.B0163:
	%7 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1434
	%8 = load i32, i32* %.ndi0002.addr, align 4, !tbaa !1401, !dbg !1434
	%9 = sext i32  %8 to i64, !dbg !1434
	%10 = load double*, double** %xa.addr, align 8, !tbaa !1398, !dbg !1434
	%11 = getelementptr double, double*  %10, i64  %9, !dbg !1434
	%12 = load double, double*  %11, align 8, !tbaa !1399, !dbg !1434
	%13 = fsub double  %7,  %12, !dbg !1434
	%14 = call double @llvm.fabs.f64 (double  %13), !dbg !1434
	call void @llvm.dbg.declare (metadata double* %dift.addr, metadata !1435, metadata !1366), !dbg !1413
	store double  %14, double* %dift.addr, align 8, !tbaa !1403, !dbg !1434
	%15 = load double, double* %dif.addr, align 8, !tbaa !1403, !dbg !1434
	%16 = fcmp uge double  %14,  %15, !dbg !1434
	br i1  %16, label %L.B0009, label %L.B0170, !dbg !1434
L.B0170:
	%17 = load i32, i32* %.ndi0002.addr, align 4, !tbaa !1401, !dbg !1434
	store i32  %17, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1434
	store double  %14, double* %dif.addr, align 8, !tbaa !1403, !dbg !1434
	br label %L.B0009
L.B0009:
	%18 = load i32, i32* %.ndi0002.addr, align 4, !tbaa !1401, !dbg !1434
	%19 = sext i32  %18 to i64, !dbg !1434
	%20 = load double*, double** %ya.addr, align 8, !tbaa !1398, !dbg !1434
	%21 = getelementptr double, double*  %20, i64  %19, !dbg !1434
	%22 = load double, double*  %21, align 8, !tbaa !1399, !dbg !1434
	call void @llvm.dbg.declare (metadata [20 x double]* %D.addr, metadata !1437, metadata !1366), !dbg !1409
	%23 = bitcast [20 x double]* %D.addr to double*, !dbg !1434
	%24 = getelementptr double, double*  %23, i64  %19, !dbg !1434
	store double  %22, double*  %24, align 8, !tbaa !1399, !dbg !1434
	%25 = load i32, i32* %.ndi0002.addr, align 4, !tbaa !1401, !dbg !1434
	%26 = sext i32  %25 to i64, !dbg !1434
	call void @llvm.dbg.declare (metadata [20 x double]* %C.addr, metadata !1438, metadata !1366), !dbg !1409
	%27 = bitcast [20 x double]* %C.addr to double*, !dbg !1434
	%28 = getelementptr double, double*  %27, i64  %26, !dbg !1434
	store double  %22, double*  %28, align 8, !tbaa !1399, !dbg !1434
	%29 = load i32, i32* %.ndi0002.addr, align 4, !tbaa !1401, !dbg !1434
	%30 = add i32  %29, 1, !dbg !1434
	store i32  %30, i32* %.ndi0002.addr, align 4, !tbaa !1401, !dbg !1434
	%31 = load i32, i32* %.lcr011025.addr, align 4, !tbaa !1401, !dbg !1434
	%32 = icmp slt i32  %30,  %31, !dbg !1434
	br i1  %32, label %L.B0163, label %L.B0164, !dbg !1434
L.B0164:
	%33 = load i32, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1439
	%34 = sext i32  %33 to i64, !dbg !1439
	%35 = load double*, double** %ya.addr, align 8, !tbaa !1398, !dbg !1439
	%36 = getelementptr double, double*  %35, i64  %34, !dbg !1439
	%37 = load double, double*  %36, align 8, !tbaa !1399, !dbg !1439
	%38 = load double*, double** %y.addr, align 8, !tbaa !1398, !dbg !1439
	store double  %37, double*  %38, align 8, !tbaa !1399, !dbg !1439

	%39 = sub i32  %33, 1, !dbg !1440
	store i32  %39, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1440
	%40 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1441
	%41 = sext i32  %40 to i64, !dbg !1441
	%42 = sub i64  %41, 1, !dbg !1441
	store i64  %42, i64* %.TRP0000.addr, align 8, !tbaa !1449, !dbg !1441
	store i64 0, i64* %.ndk0004.addr, align 8, !tbaa !1451, !dbg !1442
	%43 = icmp sle i64  %42, 0, !dbg !1442
	br i1  %43, label %L.B0166, label %L.B0171, !dbg !1442
L.B0171:
	store i64  %42, i64* %.lcr051025.addr, align 8, !tbaa !1451, !dbg !1442
	br label %L.B0165
L.B0165:
	%44 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1442
	%45 = load i64, i64* %.ndk0004.addr, align 8, !tbaa !1451, !dbg !1442
	%46 = trunc i64  %45 to i32, !dbg !1442
	%47 = sub i32  %44,  %46, !dbg !1442
	%48 = sext i32  %47 to i64, !dbg !1442
	%49 = sub i64  %48, 1, !dbg !1442
	store i64  %49, i64* %.TRP0001.addr, align 8, !tbaa !1449, !dbg !1442
	store i64 0, i64* %.ndk0006.addr, align 8, !tbaa !1451, !dbg !1442
	%50 = icmp sle i64  %49, 0, !dbg !1442
	br i1  %50, label %L.B0168, label %L.B0172, !dbg !1442
L.B0172:
	store i64  %49, i64* %.lcr051026.addr, align 8, !tbaa !1451, !dbg !1443
	store i32  %46, i32* %.r2.0254.addr, align 4, !tbaa !1401, !dbg !1443
	%51 = bitcast [20 x double]* %C.addr to i8*, !dbg !1443
	store i8*  %51, i8** %.G0001.addr, align 8, !tbaa !1398, !dbg !1443
	br label %L.B0167
L.B0167:
	%52 = load i64, i64* %.ndk0006.addr, align 8, !tbaa !1451, !dbg !1442
	%53 = trunc i64  %52 to i32, !dbg !1442
	%54 = sext i32  %53 to i64, !dbg !1442
	%55 = load double*, double** %xa.addr, align 8, !tbaa !1398, !dbg !1442
	%56 = getelementptr double, double*  %55, i64  %54, !dbg !1442
	%57 = load double, double*  %56, align 8, !tbaa !1399, !dbg !1442
	%58 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1442
	%59 = fsub double  %57,  %58, !dbg !1442
	call void @llvm.dbg.declare (metadata double* %h0.addr, metadata !1444, metadata !1366), !dbg !1417
	store double  %59, double* %h0.addr, align 8, !tbaa !1403, !dbg !1442
	%60 = load i32, i32* %.r2.0254.addr, align 4, !tbaa !1401, !dbg !1442
	%61 = load i64, i64* %.ndk0006.addr, align 8, !tbaa !1451, !dbg !1442
	%62 = trunc i64  %61 to i32, !dbg !1442
	%63 = add i32  %60,  %62, !dbg !1442
	%64 = sext i32  %63 to i64, !dbg !1442
	%65 = load double*, double** %xa.addr, align 8, !tbaa !1398, !dbg !1442
	%66 = bitcast double*  %65 to i8*, !dbg !1442
	%67 = getelementptr i8, i8*  %66, i64 8, !dbg !1442
	%68 = bitcast i8*  %67 to double*, !dbg !1442
	%69 = getelementptr double, double*  %68, i64  %64, !dbg !1442
	%70 = load double, double*  %69, align 8, !tbaa !1399, !dbg !1442
	%71 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1442
	%72 = fsub double  %70,  %71, !dbg !1442
	call void @llvm.dbg.declare (metadata double* %hp.addr, metadata !1445, metadata !1366), !dbg !1417
	store double  %72, double* %hp.addr, align 8, !tbaa !1403, !dbg !1442
	%73 = load i8*, i8** %.G0001.addr, align 8, !tbaa !1398, !dbg !1442
	%74 = getelementptr i8, i8*  %73, i64 8, !dbg !1442
	%75 = bitcast i8*  %74 to double*, !dbg !1442
	%76 = load double, double*  %75, align 8, !tbaa !1399, !dbg !1442
	%77 = load i64, i64* %.ndk0006.addr, align 8, !tbaa !1451, !dbg !1442
	%78 = trunc i64  %77 to i32, !dbg !1442
	%79 = sext i32  %78 to i64, !dbg !1442
	%80 = bitcast [20 x double]* %D.addr to double*, !dbg !1442
	%81 = getelementptr double, double*  %80, i64  %79, !dbg !1442
	%82 = load double, double*  %81, align 8, !tbaa !1399, !dbg !1442
	%83 = fsub double  %76,  %82, !dbg !1442
	call void @llvm.dbg.declare (metadata double* %W.addr, metadata !1446, metadata !1366), !dbg !1417
	store double  %83, double* %W.addr, align 8, !tbaa !1403, !dbg !1442
	%84 = load double, double* %h0.addr, align 8, !tbaa !1403, !dbg !1442
	%85 = load double, double* %hp.addr, align 8, !tbaa !1403, !dbg !1442
	%86 = fsub double  %84,  %85, !dbg !1442
	call void @llvm.dbg.declare (metadata double* %den.addr, metadata !1447, metadata !1366), !dbg !1417
	%87 = load double, double* %W.addr, align 8, !tbaa !1403, !dbg !1442
	%88 = fdiv double  %87,  %86, !dbg !1442
	%89 = load double, double* %hp.addr, align 8, !tbaa !1403, !dbg !1442
	%90 = fmul double  %89,  %88, !dbg !1442
	%91 = load i64, i64* %.ndk0006.addr, align 8, !tbaa !1451, !dbg !1442
	%92 = trunc i64  %91 to i32, !dbg !1442
	%93 = sext i32  %92 to i64, !dbg !1442
	%94 = getelementptr double, double*  %80, i64  %93, !dbg !1442
	store double  %90, double*  %94, align 8, !tbaa !1399, !dbg !1442
	%95 = load double, double* %h0.addr, align 8, !tbaa !1403, !dbg !1442
	%96 = fmul double  %95,  %88, !dbg !1442
	%97 = load i8*, i8** %.G0001.addr, align 8, !tbaa !1398, !dbg !1442
	%98 = bitcast i8*  %97 to double*, !dbg !1442
	store double  %96, double*  %98, align 8, !tbaa !1399, !dbg !1442
	%99 = add i64  %91, 1, !dbg !1442
	store i64  %99, i64* %.ndk0006.addr, align 8, !tbaa !1451, !dbg !1442
	%100 = getelementptr i8, i8*  %97, i64 8, !dbg !1443
	store i8*  %100, i8** %.G0001.addr, align 8, !tbaa !1398, !dbg !1443
	%101 = load i64, i64* %.lcr051026.addr, align 8, !tbaa !1451, !dbg !1442
	%102 = icmp slt i64  %99,  %101, !dbg !1442
	br i1  %102, label %L.B0167, label %L.B0168, !dbg !1442
L.B0168:
	%103 = load i32, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1442
	%104 = add i32  %103, 1, !dbg !1442
	%105 = mul i32  %104, 2, !dbg !1442
	%106 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1442
	%107 = load i64, i64* %.ndk0004.addr, align 8, !tbaa !1451, !dbg !1442
	%108 = trunc i64  %107 to i32, !dbg !1442
	%109 = sub i32  %106,  %108, !dbg !1442
	%110 = sub i32  %109, 1, !dbg !1442
	%111 = icmp sge i32  %105,  %110, !dbg !1442
	br i1  %111, label %L.B0018, label %L.B0173, !dbg !1442
L.B0173:
	%112 = sext i32  %103 to i64, !dbg !1442
	%113 = bitcast [20 x double]* %C.addr to i8*, !dbg !1442
	%114 = getelementptr i8, i8*  %113, i64 8, !dbg !1442
	%115 = bitcast i8*  %114 to double*, !dbg !1442
	%116 = getelementptr double, double*  %115, i64  %112, !dbg !1442
	%117 = load double, double*  %116, align 8, !tbaa !1399, !dbg !1442
	%118 = load double*, double** %dy.addr, align 8, !tbaa !1398, !dbg !1442
	store double  %117, double*  %118, align 8, !tbaa !1399, !dbg !1442
	br label %L.B0019, !dbg !1442
L.B0018:
	%119 = load i32, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1442
	%120 = sext i32  %119 to i64, !dbg !1442
	%121 = bitcast [20 x double]* %D.addr to double*, !dbg !1442
	%122 = getelementptr double, double*  %121, i64  %120, !dbg !1442
	%123 = load double, double*  %122, align 8, !tbaa !1399, !dbg !1442
	%124 = load double*, double** %dy.addr, align 8, !tbaa !1398, !dbg !1442
	store double  %123, double*  %124, align 8, !tbaa !1399, !dbg !1442

	%125 = sub i32  %119, 1, !dbg !1442
	store i32  %125, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1442
	br label %L.B0019
L.B0019:
	%126 = load double*, double** %dy.addr, align 8, !tbaa !1398, !dbg !1442
	%127 = load double, double*  %126, align 8, !tbaa !1399, !dbg !1442
	%128 = load double*, double** %y.addr, align 8, !tbaa !1398, !dbg !1442
	%129 = load double, double*  %128, align 8, !tbaa !1399, !dbg !1442
	%130 = fadd double  %127,  %129, !dbg !1442
	store double  %130, double*  %128, align 8, !tbaa !1399, !dbg !1442
	%131 = load i64, i64* %.ndk0004.addr, align 8, !tbaa !1451, !dbg !1442
	%132 = add i64  %131, 1, !dbg !1442
	store i64  %132, i64* %.ndk0004.addr, align 8, !tbaa !1451, !dbg !1442
	%133 = load i64, i64* %.lcr051025.addr, align 8, !tbaa !1451, !dbg !1442
	%134 = icmp slt i64  %132,  %133, !dbg !1442
	br i1  %134, label %L.B0165, label %L.B0166, !dbg !1442
L.B0166:
	ret void, !dbg !1441
}

%struct._ZNSt8ios_base14_Callback_listE = type <{ %struct._ZNSt8ios_base14_Callback_listE*, void (i32, %struct._ZSt8ios_base*, i32)*, i32, i32}> 
%struct._ZSt7num_putIcSt19ostreambuf_iteratorIcSt11char_traitsIcEEE = type <{ %struct.__SO__NSt6locale5facetE, [4 x i8]}> 
%struct._ZNSt6locale5facetE = type <{ i32 (...)* (...)*, i32, [4 x i8]}> 
%struct.__SO__NSt6locale5facetE = type <{ i32 (...)* (...)*, i32}> 
%struct._ZNSt8ios_base6_WordsE = type <{ i8*, i64}> 
%struct._ZSo = type <{ i32 (...)* (...)*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE}> 
%struct.__locale_data = type opaque
%struct._ZSt8ios_base = type <{ i32 (...)* (...)*, i64, i64, i32, i32, i32, [4 x i8], %struct._ZNSt8ios_base14_Callback_listE*, %struct._ZNSt8ios_base6_WordsE, [8 x %struct._ZNSt8ios_base6_WordsE], i32, [4 x i8], %struct._ZNSt8ios_base6_WordsE*, %struct._ZSt6locale}> 
%struct._ZSt5ctypeIcE = type <{ %struct.__SO__NSt6locale5facetE, [4 x i8], %struct.__locale_struct*, i8, [7 x i8], i32*, i32*, i16*, i8, [256 x i8], [256 x i8], i8, [6 x i8]}> 
%struct._ZSt9basic_iosIcSt11char_traitsIcEE = type <{ %struct._ZSt8ios_base, %struct._ZSo*, i8, i8, [6 x i8], %struct._ZSt15basic_streambufIcSt11char_traitsIcEE*, %struct._ZSt5ctypeIcE*, %struct._ZSt7num_putIcSt19ostreambuf_iteratorIcSt11char_traitsIcEEE*, %struct._ZSt7num_getIcSt19istreambuf_iteratorIcSt11char_traitsIcEEE*}> 
%struct._ZSt6locale = type <{ %struct._ZNSt6locale5_ImplE*}> 
%struct._ZSt15basic_streambufIcSt11char_traitsIcEE = type <{ i32 (...)* (...)*, i8*, i8*, i8*, i8*, i8*, i8*, %struct._ZSt6locale}> 
%struct.__locale_struct = type <{ [13 x %struct.__locale_data*], i16*, i32*, i32*, [13 x i8*]}> 
%struct._ZSt7num_getIcSt19istreambuf_iteratorIcSt11char_traitsIcEEE = type <{ %struct.__SO__NSt6locale5facetE, [4 x i8]}> 
%struct._ZNSt6locale5_ImplE = type <{ i32, [4 x i8], %struct._ZNSt6locale5facetE**, i64, %struct._ZNSt6locale5facetE**, i8**}> 

define internal void @_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a93656ratintEPKdS1_idRdS2_(double* %xa.arg, double* %ya.arg, i32 signext %n.arg, double %x.arg, double* %y.arg, double* %dy.arg) #0 !dbg !1453 {
L.entry:
	%xa.addr = alloca double*, align 8
	%ya.addr = alloca double*, align 8
	%n.addr = alloca i32, align 4
	%x.addr = alloca double, align 8
	%y.addr = alloca double*, align 8
	%dy.addr = alloca double*, align 8
	%ns.addr = alloca i32, align 4
	%hh.addr = alloca double, align 8
	%i.addr = alloca i32, align 4
	%.lcr012049.addr = alloca i32, align 4
	%h.addr = alloca double, align 8
	%C.addr = alloca [20 x double], align 8
	%D.addr = alloca [20 x double], align 8
	%.D0007.addr = alloca i32, align 4
	%m1.addr = alloca i32, align 4
	%.x2049.addr = alloca i32, align 4
	%.lcr012050.addr = alloca i32, align 4
	%.G0011.addr = alloca i8*, align 8
	%w.addr = alloca double, align 8
	%t.addr = alloca double, align 8
	%dd.addr = alloca double, align 8
	%..inline.addr = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.1 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.2 = alloca i32, align 4
	%..inline.addr.3 = alloca i64, align 8
	%.Q0003.addr = alloca %struct._ZSo*, align 8
	%.I0000.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata double** %xa.addr, metadata !1477, metadata !1366), !dbg !1454
	store double* %xa.arg, double** %xa.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %xa.addr, metadata !1478, metadata !1366), !dbg !1454
	call void @llvm.dbg.declare (metadata double** %ya.addr, metadata !1479, metadata !1366), !dbg !1454
	store double* %ya.arg, double** %ya.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %ya.addr, metadata !1480, metadata !1366), !dbg !1454
	call void @llvm.dbg.declare (metadata i32* %n.addr, metadata !1481, metadata !1366), !dbg !1454
	store i32 %n.arg, i32* %n.addr, align 4, !tbaa !1399
	call void @llvm.dbg.declare (metadata i32* %n.addr, metadata !1482, metadata !1366), !dbg !1454
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !1483, metadata !1366), !dbg !1454
	store double %x.arg, double* %x.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !1484, metadata !1366), !dbg !1454
	call void @llvm.dbg.declare (metadata double** %y.addr, metadata !1485, metadata !1366), !dbg !1454
	store double* %y.arg, double** %y.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %y.addr, metadata !1486, metadata !1366), !dbg !1454
	call void @llvm.dbg.declare (metadata double** %dy.addr, metadata !1487, metadata !1366), !dbg !1454
	store double* %dy.arg, double** %dy.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata double** %dy.addr, metadata !1488, metadata !1366), !dbg !1454
	call void @llvm.dbg.declare (metadata i32* %ns.addr, metadata !1490, metadata !1366), !dbg !1454
	store i32 0, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1489
	%0 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1491
	%1 = load double*, double** %xa.addr, align 8, !tbaa !1398, !dbg !1491
	%2 = load double, double*  %1, align 8, !tbaa !1399, !dbg !1491
	%3 = fsub double  %0,  %2, !dbg !1491
	%4 = call double @llvm.fabs.f64 (double  %3), !dbg !1491
	call void @llvm.dbg.declare (metadata double* %hh.addr, metadata !1492, metadata !1366), !dbg !1454
	store double  %4, double* %hh.addr, align 8, !tbaa !1403, !dbg !1491
	call void @llvm.dbg.declare (metadata i32* %i.addr, metadata !1494, metadata !1366), !dbg !1454
	store i32 0, i32* %i.addr, align 4, !tbaa !1401, !dbg !1493
	%5 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1493
	%6 = icmp sle i32  %5, 0, !dbg !1493
	br i1  %6, label %L.B0022, label %L.B0176, !dbg !1493
L.B0176:
	store i32  %5, i32* %.lcr012049.addr, align 4, !tbaa !1401, !dbg !1493
	br label %L.B0021
L.B0021:
	%7 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1495
	%8 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1495
	%9 = sext i32  %8 to i64, !dbg !1495
	%10 = load double*, double** %xa.addr, align 8, !tbaa !1398, !dbg !1495
	%11 = getelementptr double, double*  %10, i64  %9, !dbg !1495
	%12 = load double, double*  %11, align 8, !tbaa !1399, !dbg !1495
	%13 = fsub double  %7,  %12, !dbg !1495
	%14 = call double @llvm.fabs.f64 (double  %13), !dbg !1495
	call void @llvm.dbg.declare (metadata double* %h.addr, metadata !1496, metadata !1366), !dbg !1454
	store double  %14, double* %h.addr, align 8, !tbaa !1403, !dbg !1495
	%15 = fcmp une double  %14,  0.00000000000000000E+0, !dbg !1497
	br i1  %15, label %L.B0023, label %L.B0177, !dbg !1497
L.B0177:
	%16 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1498
	%17 = sext i32  %16 to i64, !dbg !1498
	%18 = load double*, double** %ya.addr, align 8, !tbaa !1398, !dbg !1498
	%19 = getelementptr double, double*  %18, i64  %17, !dbg !1498
	%20 = load double, double*  %19, align 8, !tbaa !1399, !dbg !1498
	%21 = load double*, double** %y.addr, align 8, !tbaa !1398, !dbg !1498
	store double  %20, double*  %21, align 8, !tbaa !1399, !dbg !1498
	%22 = load double*, double** %dy.addr, align 8, !tbaa !1398, !dbg !1499
	store double  0.00000000000000000E+0, double*  %22, align 8, !tbaa !1399, !dbg !1499
	br label %L.R0002, !dbg !1500
L.B0023:
	%23 = load double, double* %h.addr, align 8, !tbaa !1403, !dbg !1501
	%24 = load double, double* %hh.addr, align 8, !tbaa !1403, !dbg !1501
	%25 = fcmp uge double  %23,  %24, !dbg !1501
	br i1  %25, label %L.B0025, label %L.B0178, !dbg !1501
L.B0178:
	%26 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1502
	store i32  %26, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1502
	store double  %23, double* %hh.addr, align 8, !tbaa !1403, !dbg !1503
	br label %L.B0025
L.B0025:
	%27 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1504
	%28 = sext i32  %27 to i64, !dbg !1504
	%29 = load double*, double** %ya.addr, align 8, !tbaa !1398, !dbg !1504
	%30 = getelementptr double, double*  %29, i64  %28, !dbg !1504
	%31 = load double, double*  %30, align 8, !tbaa !1399, !dbg !1504
	call void @llvm.dbg.declare (metadata [20 x double]* %C.addr, metadata !1505, metadata !1366), !dbg !1454
	%32 = bitcast [20 x double]* %C.addr to double*, !dbg !1504
	%33 = getelementptr double, double*  %32, i64  %28, !dbg !1504
	store double  %31, double*  %33, align 8, !tbaa !1399, !dbg !1504
	%34 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1506
	%35 = sext i32  %34 to i64, !dbg !1506
	%36 = load double*, double** %ya.addr, align 8, !tbaa !1398, !dbg !1506
	%37 = getelementptr double, double*  %36, i64  %35, !dbg !1506
	%38 = load double, double*  %37, align 8, !tbaa !1399, !dbg !1506
	%39 = fadd double  %38,  1.00000000000000008E-30, !dbg !1506
	call void @llvm.dbg.declare (metadata [20 x double]* %D.addr, metadata !1507, metadata !1366), !dbg !1454
	%40 = bitcast [20 x double]* %D.addr to double*, !dbg !1506
	%41 = getelementptr double, double*  %40, i64  %35, !dbg !1506
	store double  %39, double*  %41, align 8, !tbaa !1399, !dbg !1506
	%42 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1508

	%43 = add i32  %42, 1, !dbg !1508
	store i32  %43, i32* %i.addr, align 4, !tbaa !1401, !dbg !1508
	%44 = load i32, i32* %.lcr012049.addr, align 4, !tbaa !1401, !dbg !1508
	%45 = icmp slt i32  %43,  %44, !dbg !1508
	br i1  %45, label %L.B0021, label %L.B0022, !dbg !1508
L.B0022:
	%46 = load i32, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1509

	%47 = sub i32  %46, 1, !dbg !1509
	store i32  %47, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1509
	%48 = sext i32  %46 to i64, !dbg !1509
	%49 = load double*, double** %ya.addr, align 8, !tbaa !1398, !dbg !1509
	%50 = getelementptr double, double*  %49, i64  %48, !dbg !1509
	%51 = load double, double*  %50, align 8, !tbaa !1399, !dbg !1509
	%52 = load double*, double** %y.addr, align 8, !tbaa !1398, !dbg !1509
	store double  %51, double*  %52, align 8, !tbaa !1399, !dbg !1509
	call void @llvm.dbg.declare (metadata i32* %m1.addr, metadata !1511, metadata !1366), !dbg !1454
	store i32 1, i32* %m1.addr, align 4, !tbaa !1401, !dbg !1510
	%53 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1510
	%54 = icmp sge i32 1,  %53, !dbg !1510
	br i1  %54, label %L.B0027, label %L.B0179, !dbg !1510
L.B0179:
	%55 = sub i32  %53, 1, !dbg !1510
	store i32  %55, i32* %.x2049.addr, align 4, !tbaa !1401, !dbg !1510
	br label %L.B0026
L.B0026:
	store i32 0, i32* %i.addr, align 4, !tbaa !1401, !dbg !1512
	%56 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1512
	%57 = load i32, i32* %m1.addr, align 4, !tbaa !1401, !dbg !1512
	%58 = sub i32  %56,  %57, !dbg !1512
	%59 = icmp sle i32  %58, 0, !dbg !1512
	br i1  %59, label %L.B0029, label %L.B0180, !dbg !1512
L.B0180:
	store i32  %58, i32* %.lcr012050.addr, align 4, !tbaa !1401, !dbg !1512
	%60 = bitcast [20 x double]* %C.addr to i8*, !dbg !1512
	%61 = getelementptr i8, i8*  %60, i64 8, !dbg !1512
	store i8*  %61, i8** %.G0011.addr, align 8, !tbaa !1398, !dbg !1512
	br label %L.B0028
L.B0028:
	%62 = load i8*, i8** %.G0011.addr, align 8, !tbaa !1398, !dbg !1513
	%63 = bitcast i8*  %62 to double*, !dbg !1513
	%64 = load double, double*  %63, align 8, !tbaa !1399, !dbg !1513
	%65 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1513
	%66 = sext i32  %65 to i64, !dbg !1513
	%67 = bitcast [20 x double]* %D.addr to double*, !dbg !1513
	%68 = getelementptr double, double*  %67, i64  %66, !dbg !1513
	%69 = load double, double*  %68, align 8, !tbaa !1399, !dbg !1513
	%70 = fsub double  %64,  %69, !dbg !1513
	call void @llvm.dbg.declare (metadata double* %w.addr, metadata !1514, metadata !1366), !dbg !1454
	store double  %70, double* %w.addr, align 8, !tbaa !1403, !dbg !1513
	%71 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1515
	%72 = load i32, i32* %m1.addr, align 4, !tbaa !1401, !dbg !1515
	%73 = add i32  %71,  %72, !dbg !1515
	%74 = sext i32  %73 to i64, !dbg !1515
	%75 = load double*, double** %xa.addr, align 8, !tbaa !1398, !dbg !1515
	%76 = getelementptr double, double*  %75, i64  %74, !dbg !1515
	%77 = load double, double*  %76, align 8, !tbaa !1399, !dbg !1515
	%78 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !1515
	%79 = fsub double  %77,  %78, !dbg !1515
	%80 = sext i32  %71 to i64, !dbg !1516
	%81 = getelementptr double, double*  %67, i64  %80, !dbg !1516
	%82 = load double, double*  %81, align 8, !tbaa !1399, !dbg !1516
	%83 = getelementptr double, double*  %75, i64  %80, !dbg !1516
	%84 = load double, double*  %83, align 8, !tbaa !1399, !dbg !1516
	%85 = fsub double  %84,  %78, !dbg !1516
	%86 = fmul double  %82,  %85, !dbg !1516
	%87 = fdiv double  %86,  %79, !dbg !1516
	call void @llvm.dbg.declare (metadata double* %t.addr, metadata !1517, metadata !1366), !dbg !1454
	store double  %87, double* %t.addr, align 8, !tbaa !1403, !dbg !1516
	%88 = load i8*, i8** %.G0011.addr, align 8, !tbaa !1398, !dbg !1518
	%89 = bitcast i8*  %88 to double*, !dbg !1518
	%90 = load double, double*  %89, align 8, !tbaa !1399, !dbg !1518
	%91 = fsub double  %87,  %90, !dbg !1518
	call void @llvm.dbg.declare (metadata double* %dd.addr, metadata !1519, metadata !1366), !dbg !1454
	store double  %91, double* %dd.addr, align 8, !tbaa !1403, !dbg !1518
	%92 = fcmp une double  %91,  0.00000000000000000E+0, !dbg !1520
	br i1  %92, label %L.B0030, label %L.B0181, !dbg !1520
L.B0181:
	%93 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !1521
	%94 = icmp ne i8*  %93,  null, !dbg !1521
	br i1  %94, label %L..inline.10044, label %L.B0182, !dbg !1521
L.B0182:
	%95 = bitcast %struct._ZSo* @_ZSt4cerr to i8*, !dbg !1522
	%96 = bitcast %struct._ZSo* @_ZSt4cerr to i8**, !dbg !1522
	%97 = load i8*, i8**  %96, align 8, !tbaa !1398, !dbg !1522
	%98 = getelementptr i8, i8*  %97, i64 18446744073709551592, !dbg !1522
	%99 = bitcast i8*  %98 to i64*, !dbg !1522
	%100 = load i64, i64*  %99, align 8, !tbaa !1399, !dbg !1522
	%101 = getelementptr i8, i8*  %95, i64  %100, !dbg !1522
	%102 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr to i8**, !dbg !1522
	store i8*  %101, i8**  %102, align 8, !tbaa !1398, !dbg !1522
	%103 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr, align 8, !tbaa !1398, !dbg !1522
	%104 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %103 to i8*, !dbg !1522
	%105 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.1 to i8**, !dbg !1522
	store i8*  %104, i8**  %105, align 8, !tbaa !1398, !dbg !1522
	%106 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.1, align 8, !tbaa !1398, !dbg !1522
	%107 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %106 to i8*, !dbg !1522
	%108 = getelementptr i8, i8*  %107, i64 32, !dbg !1522
	%109 = bitcast i8*  %108 to i32*, !dbg !1522
	%110 = load i32, i32*  %109, align 4, !tbaa !1399, !dbg !1522
	%111 = or i32  %110, 1, !dbg !1522
	call void  @_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate (%struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %103, i32  %111), !dbg !1522
	br label %L..inline.10065, !dbg !1522
L..inline.10044:
	%112 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !1525
	%113 = call i64  @strlen (i8*  %112) nounwind, !dbg !1525
	%114 = call %struct._ZSo*  @_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l (%struct._ZSo* @_ZSt4cerr, i8*  %112, i64  %113), !dbg !1525
	%115 = bitcast %struct._ZSo*  %114 to i8*, !dbg !1525
	%116 = bitcast %struct._ZSo** %.Q0003.addr to i8**, !dbg !1525
	store i8*  %115, i8**  %116, align 8, !tbaa !1398, !dbg !1525
	br label %L..inline.10065
L..inline.10065:
	call void  @exit (i32 111) nounwind noreturn, !dbg !1526
	br label %L.B0030
L.B0030:
	%117 = load double, double* %w.addr, align 8, !tbaa !1403, !dbg !1527
	%118 = load double, double* %dd.addr, align 8, !tbaa !1403, !dbg !1527
	%119 = fdiv double  %117,  %118, !dbg !1527
	%120 = load i8*, i8** %.G0011.addr, align 8, !tbaa !1398, !dbg !1528
	%121 = bitcast i8*  %120 to double*, !dbg !1528
	%122 = load double, double*  %121, align 8, !tbaa !1399, !dbg !1528
	%123 = fmul double  %119,  %122, !dbg !1528
	%124 = load i32, i32* %i.addr, align 4, !tbaa !1401, !dbg !1528
	%125 = sext i32  %124 to i64, !dbg !1528
	%126 = bitcast [20 x double]* %D.addr to double*, !dbg !1528
	%127 = getelementptr double, double*  %126, i64  %125, !dbg !1528
	store double  %123, double*  %127, align 8, !tbaa !1399, !dbg !1528
	%128 = load double, double* %t.addr, align 8, !tbaa !1403, !dbg !1529
	%129 = fmul double  %128,  %119, !dbg !1529
	%130 = getelementptr i8, i8*  %120, i64 18446744073709551608, !dbg !1529
	%131 = bitcast i8*  %130 to double*, !dbg !1529
	store double  %129, double*  %131, align 8, !tbaa !1399, !dbg !1529

	%132 = add i32  %124, 1, !dbg !1530
	store i32  %132, i32* %i.addr, align 4, !tbaa !1401, !dbg !1530
	%133 = getelementptr i8, i8*  %120, i64 8, !dbg !1512
	store i8*  %133, i8** %.G0011.addr, align 8, !tbaa !1398, !dbg !1512
	%134 = load i32, i32* %.lcr012050.addr, align 4, !tbaa !1401, !dbg !1530
	%135 = icmp slt i32  %132,  %134, !dbg !1530
	br i1  %135, label %L.B0028, label %L.B0029, !dbg !1530
L.B0029:
	%136 = load i32, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1531
	%137 = add i32  %136, 1, !dbg !1531
	%138 = mul i32  %137, 2, !dbg !1531
	%139 = load i32, i32* %n.addr, align 4, !tbaa !1401, !dbg !1531
	%140 = load i32, i32* %m1.addr, align 4, !tbaa !1401, !dbg !1531
	%141 = sub i32  %139,  %140, !dbg !1531
	%142 = icmp sge i32  %138,  %141, !dbg !1531
	br i1  %142, label %L.B0031, label %L.B0183, !dbg !1531
L.B0183:
	%143 = sext i32  %136 to i64, !dbg !1531
	%144 = bitcast [20 x double]* %C.addr to i8*, !dbg !1531
	%145 = getelementptr i8, i8*  %144, i64 8, !dbg !1531
	%146 = bitcast i8*  %145 to double*, !dbg !1531
	%147 = getelementptr double, double*  %146, i64  %143, !dbg !1531
	%148 = load double, double*  %147, align 8, !tbaa !1399, !dbg !1531
	store double  %148, double* %.I0000.addr, align 8, !tbaa !1403, !dbg !1531
	br label %L.B0032, !dbg !1531
L.B0031:
	%149 = load i32, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1531

	%150 = sub i32  %149, 1, !dbg !1531
	store i32  %150, i32* %ns.addr, align 4, !tbaa !1401, !dbg !1531
	%151 = sext i32  %149 to i64, !dbg !1531
	%152 = bitcast [20 x double]* %D.addr to double*, !dbg !1531
	%153 = getelementptr double, double*  %152, i64  %151, !dbg !1531
	%154 = load double, double*  %153, align 8, !tbaa !1399, !dbg !1531
	store double  %154, double* %.I0000.addr, align 8, !tbaa !1403, !dbg !1531
	br label %L.B0032
L.B0032:
	%155 = load double, double* %.I0000.addr, align 8, !tbaa !1403, !dbg !1531
	%156 = load double*, double** %dy.addr, align 8, !tbaa !1398, !dbg !1531
	store double  %155, double*  %156, align 8, !tbaa !1399, !dbg !1531
	%157 = load double*, double** %y.addr, align 8, !tbaa !1398, !dbg !1531
	%158 = load double, double*  %157, align 8, !tbaa !1399, !dbg !1531
	%159 = fadd double  %158,  %155, !dbg !1531
	store double  %159, double*  %157, align 8, !tbaa !1399, !dbg !1531
	%160 = load i32, i32* %m1.addr, align 4, !tbaa !1401, !dbg !1532

	%161 = add i32  %160, 1, !dbg !1532
	store i32  %161, i32* %m1.addr, align 4, !tbaa !1401, !dbg !1532
	%162 = load i32, i32* %.x2049.addr, align 4, !tbaa !1401, !dbg !1510
	%163 = sub i32  %162, 1, !dbg !1510
	store i32  %163, i32* %.x2049.addr, align 4, !tbaa !1401, !dbg !1510
	%164 = icmp sgt i32  %163, 0, !dbg !1532
	br i1  %164, label %L.B0026, label %L.B0027, !dbg !1532
L.B0027:
	br label %L.N0002
L.N0002:
	br label %L.R0002
L.R0002:
	ret void, !dbg !1533
}
define double @_Z14Romberg_simpleRK11T1DFunctionddd(%struct.T1DFunction* %func.arg, double %a.arg, double %b.arg, double %absacc.arg) #0 !dbg !1537 {
L.entry:
	%func.addr = alloca %struct.T1DFunction*, align 8
	%a.addr = alloca double, align 8
	%b.addr = alloca double, align 8
	%absacc.addr = alloca double, align 8
	%it.addr = alloca i32, align 4
	%H.addr = alloca [9 x double], align 8
	%dresult.addr = alloca double, align 8
	%j.addr = alloca i32, align 4
	%.x3073.addr = alloca i32, align 4
	%S.addr = alloca [9 x double], align 8
	%.G0035.addr = alloca i8*, align 8
	%.G0033.addr = alloca i8*, align 8
	%..inline.addr = alloca %struct.T1DFunction*, align 8
	%..inline.addr.1 = alloca double, align 8
	%..inline.addr.2 = alloca double, align 8
	%..inline.addr.3 = alloca double*, align 8
	%.Q0004.addr = alloca double, align 8
	%.Q0005.addr = alloca double, align 8
	%..inline.addr.4 = alloca double, align 8
	%..inline.addr.5 = alloca double, align 8
	%..inline.addr.6 = alloca double, align 8
	%.x3074.addr = alloca i32, align 4
	%.Q0006.addr = alloca double, align 8
	%result.addr = alloca double, align 8
	%k1.addr = alloca i32, align 4
	%..inline.addr.7 = alloca double*, align 8
	%..inline.addr.8 = alloca double*, align 8
	%..inline.addr.9 = alloca i32, align 4
	%..inline.addr.320 = alloca i32, align 4
	%..inline.addr.330 = alloca double, align 8
	%.ndi0006.addr = alloca i32, align 4
	%.lcr013074.addr = alloca i32, align 4
	%..inline.addr.360 = alloca double, align 8
	%..inline.addr.370 = alloca [20 x double], align 8
	%..inline.addr.380 = alloca [20 x double], align 8
	%.TRP0002.addr = alloca i64, align 8
	%.ndk0012.addr = alloca i64, align 8
	%.lcr053074.addr = alloca i64, align 8
	%.TRP0003.addr = alloca i64, align 8
	%.ndk0014.addr = alloca i64, align 8
	%.lcr053075.addr = alloca i64, align 8
	%.r4.0386.addr = alloca i32, align 4
	%.G0020.addr = alloca i8*, align 8
	%..inline.addr.470 = alloca double, align 8
	%..inline.addr.480 = alloca double, align 8
	%..inline.addr.490 = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %func.addr, metadata !1565, metadata !1366), !dbg !1538
	store %struct.T1DFunction* %func.arg, %struct.T1DFunction** %func.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %func.addr, metadata !1566, metadata !1366), !dbg !1538
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1567, metadata !1366), !dbg !1538
	store double %a.arg, double* %a.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1568, metadata !1366), !dbg !1538
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1569, metadata !1366), !dbg !1538
	store double %b.arg, double* %b.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1570, metadata !1366), !dbg !1538
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !1571, metadata !1366), !dbg !1538
	store double %absacc.arg, double* %absacc.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !1572, metadata !1366), !dbg !1538
	call void @llvm.dbg.declare (metadata i32* %it.addr, metadata !1574, metadata !1366), !dbg !1538
	store i32 0, i32* %it.addr, align 4, !tbaa !1401, !dbg !1573
	call void @llvm.dbg.declare (metadata [9 x double]* %H.addr, metadata !1579, metadata !1366), !dbg !1538
	%0 = bitcast [9 x double]* %H.addr to double*, !dbg !1575
	store double  1.00000000000000000E+0, double*  %0, align 8, !tbaa !1399, !dbg !1575
	call void @llvm.dbg.declare (metadata double* %dresult.addr, metadata !1581, metadata !1366), !dbg !1538
	store double  0.00000000000000000E+0, double* %dresult.addr, align 8, !tbaa !1403, !dbg !1580
	call void @llvm.dbg.declare (metadata i32* %j.addr, metadata !1583, metadata !1366), !dbg !1538
	store i32 1, i32* %j.addr, align 4, !tbaa !1401, !dbg !1582
	store i32 8, i32* %.x3073.addr, align 4, !tbaa !1401, !dbg !1582
	call void @llvm.dbg.declare (metadata [9 x double]* %S.addr, metadata !1584, metadata !1366), !dbg !1538
	%1 = bitcast [9 x double]* %S.addr to i8*, !dbg !1582
	%2 = getelementptr i8, i8*  %1, i64 8, !dbg !1582
	store i8*  %2, i8** %.G0035.addr, align 8, !tbaa !1398, !dbg !1582
	%3 = bitcast [9 x double]* %H.addr to i8*, !dbg !1582
	%4 = getelementptr i8, i8*  %3, i64 8, !dbg !1582
	store i8*  %4, i8** %.G0033.addr, align 8, !tbaa !1398, !dbg !1582
	br label %L.B0033
L.B0033:
	%5 = load %struct.T1DFunction*, %struct.T1DFunction** %func.addr, align 8, !tbaa !1398, !dbg !1585
	%6 = bitcast %struct.T1DFunction*  %5 to i8*, !dbg !1585
	%7 = bitcast %struct.T1DFunction** %..inline.addr to i8**, !dbg !1585
	store i8*  %6, i8**  %7, align 8, !tbaa !1398, !dbg !1585
	%8 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !1585
	store double  %8, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1585
	%9 = load double, double* %b.addr, align 8, !tbaa !1403, !dbg !1585
	store double  %9, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !1585
	%10 = bitcast [9 x double]* %S.addr to i8*, !dbg !1585
	%11 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1585
	%12 = sext i32  %11 to i64, !dbg !1585
	%13 = sub i64  %12, 1, !dbg !1585
	%14 = mul i64  %13, 8, !dbg !1585
	%15 = getelementptr i8, i8*  %10, i64  %14, !dbg !1585
	%16 = bitcast double** %..inline.addr.3 to i8**, !dbg !1585
	store i8*  %15, i8**  %16, align 8, !tbaa !1398, !dbg !1585
	%17 = icmp ne i32  %11, 1, !dbg !1585
	br i1  %17, label %L..inline.10296, label %L.B0203, !dbg !1585
L.B0203:
	%18 = load %struct.T1DFunction*, %struct.T1DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !1586
	%19 = bitcast %struct.T1DFunction*  %18 to double (%struct.T1DFunction*, double)***, !dbg !1586
	%20 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %19, align 8, !tbaa !1398, !dbg !1586
	%21 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %20, align 8, !tbaa !1398, !dbg !1586
	%22 = call double  %21 (%struct.T1DFunction*  %18, double  %9), !dbg !1586
	store double  %22, double* %.Q0004.addr, align 8, !tbaa !1403, !dbg !1586
	%23 = load %struct.T1DFunction*, %struct.T1DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !1586
	%24 = load double, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1586
	%25 = bitcast %struct.T1DFunction*  %23 to double (%struct.T1DFunction*, double)***, !dbg !1586
	%26 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %25, align 8, !tbaa !1398, !dbg !1586
	%27 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %26, align 8, !tbaa !1398, !dbg !1586
	%28 = call double  %27 (%struct.T1DFunction*  %23, double  %24), !dbg !1586
	store double  %28, double* %.Q0005.addr, align 8, !tbaa !1403, !dbg !1586
	%29 = load double, double* %.Q0004.addr, align 8, !tbaa !1403, !dbg !1586
	%30 = fadd double  %29,  %28, !dbg !1586
	%31 = load double, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !1586
	%32 = load double, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1586
	%33 = fsub double  %31,  %32, !dbg !1586
	%34 = fmul double  %33,  5.00000000000000000E-1, !dbg !1586
	%35 = fmul double  %30,  %34, !dbg !1586
	%36 = load double*, double** %..inline.addr.3, align 8, !tbaa !1398, !dbg !1586
	store double  %35, double*  %36, align 8, !tbaa !1399, !dbg !1586
	store i32 1, i32* %it.addr, align 4, !tbaa !1401, !dbg !1586
	br label %L..inline.10302, !dbg !1587
L..inline.10296:
	%37 = load double, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !1588
	%38 = load double, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1588
	%39 = fsub double  %37,  %38, !dbg !1588
	%40 = load i32, i32* %it.addr, align 4, !tbaa !1401, !dbg !1588
	%41 = sitofp i32  %40 to double, !dbg !1588
	%42 = fdiv double  %39,  %41, !dbg !1588
	store double  %42, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !1588
	%43 = call double @llvm.fma.f64 (double  %42, double  5.00000000000000000E-1, double  %38), !dbg !1589
	store double  %43, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1589
	store double  0.00000000000000000E+0, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1590
	%44 = icmp sle i32  %40, 0, !dbg !1591
	br i1  %44, label %L..inline.10311, label %L.B0204, !dbg !1591
L.B0204:
	store i32  %40, i32* %.x3074.addr, align 4, !tbaa !1401, !dbg !1591
	br label %L..inline.10310
L..inline.10310:
	%45 = load %struct.T1DFunction*, %struct.T1DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !1592
	%46 = load double, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1592
	%47 = bitcast %struct.T1DFunction*  %45 to double (%struct.T1DFunction*, double)***, !dbg !1592
	%48 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %47, align 8, !tbaa !1398, !dbg !1592
	%49 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %48, align 8, !tbaa !1398, !dbg !1592
	%50 = call double  %49 (%struct.T1DFunction*  %45, double  %46), !dbg !1592
	store double  %50, double* %.Q0006.addr, align 8, !tbaa !1403, !dbg !1592
	%51 = load double, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1592
	%52 = fadd double  %50,  %51, !dbg !1592
	store double  %52, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1592
	%53 = load double, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !1592
	%54 = load double, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1592
	%55 = fadd double  %53,  %54, !dbg !1592
	store double  %55, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1592
	%56 = load i32, i32* %.x3074.addr, align 4, !tbaa !1401, !dbg !1591
	%57 = sub i32  %56, 1, !dbg !1591
	store i32  %57, i32* %.x3074.addr, align 4, !tbaa !1401, !dbg !1591
	%58 = icmp sgt i32  %57, 0, !dbg !1593
	br i1  %58, label %L..inline.10310, label %L..inline.10311, !dbg !1593
L..inline.10311:
	%59 = load double, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !1593
	%60 = load double, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1593
	%61 = fsub double  %59,  %60, !dbg !1593
	%62 = load double, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1593
	%63 = fmul double  %61,  %62, !dbg !1593
	%64 = load i32, i32* %it.addr, align 4, !tbaa !1401, !dbg !1593
	%65 = sitofp i32  %64 to double, !dbg !1593
	%66 = fdiv double  %63,  %65, !dbg !1593
	%67 = load double*, double** %..inline.addr.3, align 8, !tbaa !1398, !dbg !1593
	%68 = load double, double*  %67, align 8, !tbaa !1399, !dbg !1593
	%69 = fadd double  %66,  %68, !dbg !1593
	%70 = fmul double  %69,  5.00000000000000000E-1, !dbg !1593
	store double  %70, double*  %67, align 8, !tbaa !1399, !dbg !1593
	%71 = mul i32  %64, 2, !dbg !1594
	store i32  %71, i32* %it.addr, align 4, !tbaa !1401, !dbg !1594
	br label %L..inline.10302
L..inline.10302:
	%72 = load i8*, i8** %.G0035.addr, align 8, !tbaa !1398, !dbg !1595
	%73 = getelementptr i8, i8*  %72, i64 18446744073709551608, !dbg !1595
	%74 = bitcast i8*  %73 to double*, !dbg !1595
	%75 = load double, double*  %74, align 8, !tbaa !1399, !dbg !1595
	call void @llvm.dbg.declare (metadata double* %result.addr, metadata !1596, metadata !1366), !dbg !1538
	store double  %75, double* %result.addr, align 8, !tbaa !1403, !dbg !1595
	%76 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1597
	%77 = icmp sle i32  %76, 0, !dbg !1597
	br i1  %77, label %L.B0035, label %L.B0205, !dbg !1597
L.B0205:
	%78 = icmp slt i32  %76, 3, !dbg !1598
	%79 = select i1  %78, i32  %76, i32 3, !dbg !1598
	call void @llvm.dbg.declare (metadata i32* %k1.addr, metadata !1599, metadata !1366), !dbg !1538
	%80 = bitcast [9 x double]* %H.addr to i8*, !dbg !1600
	%81 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1600
	%82 = sub i32  %81,  %79, !dbg !1600
	%83 = sext i32  %82 to i64, !dbg !1600
	%84 = mul i64  %83, 8, !dbg !1600
	%85 = getelementptr i8, i8*  %80, i64  %84, !dbg !1600
	%86 = bitcast double** %..inline.addr.7 to i8**, !dbg !1600
	store i8*  %85, i8**  %86, align 8, !tbaa !1398, !dbg !1600
	%87 = bitcast [9 x double]* %S.addr to i8*, !dbg !1600
	%88 = mul i64  %83, 8, !dbg !1600
	%89 = getelementptr i8, i8*  %87, i64  %88, !dbg !1600
	%90 = bitcast double** %..inline.addr.8 to i8**, !dbg !1600
	store i8*  %89, i8**  %90, align 8, !tbaa !1398, !dbg !1600
	store i32  %79, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !1600
	store i32 0, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1601
	%91 = load double*, double** %..inline.addr.7, align 8, !tbaa !1398, !dbg !1602
	%92 = load double, double*  %91, align 8, !tbaa !1399, !dbg !1602
	%93 = fsub double -0.00000000e+00,  %92, !dbg !1602
	%94 = call double @llvm.fabs.f64 (double  %93), !dbg !1602
	store double  %94, double* %..inline.addr.330, align 8, !tbaa !1403, !dbg !1602
	store i32 0, i32* %.ndi0006.addr, align 4, !tbaa !1401, !dbg !1603
	%95 = icmp slt i32  %81, 3, !dbg !1603
	%96 = select i1  %95, i32  %81, i32 3, !dbg !1603
	%97 = icmp sle i32  %96, 0, !dbg !1603
	br i1  %97, label %L.B0198, label %L.B0206, !dbg !1603
L.B0206:
	%98 = icmp slt i32  %81, 3, !dbg !1604
	%99 = select i1  %98, i32  %81, i32 3, !dbg !1604
	store i32  %99, i32* %.lcr013074.addr, align 4, !tbaa !1401, !dbg !1604
	br label %L.B0197
L.B0197:
	%100 = load i32, i32* %.ndi0006.addr, align 4, !tbaa !1401, !dbg !1603
	%101 = sext i32  %100 to i64, !dbg !1603
	%102 = load double*, double** %..inline.addr.7, align 8, !tbaa !1398, !dbg !1603
	%103 = getelementptr double, double*  %102, i64  %101, !dbg !1603
	%104 = load double, double*  %103, align 8, !tbaa !1399, !dbg !1603
	%105 = fsub double -0.00000000e+00,  %104, !dbg !1603
	%106 = call double @llvm.fabs.f64 (double  %105), !dbg !1603
	store double  %106, double* %..inline.addr.360, align 8, !tbaa !1403, !dbg !1603
	%107 = load double, double* %..inline.addr.330, align 8, !tbaa !1403, !dbg !1603
	%108 = fcmp uge double  %106,  %107, !dbg !1603
	br i1  %108, label %L..inline.10329, label %L.B0207, !dbg !1603
L.B0207:
	store i32  %100, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1603
	store double  %106, double* %..inline.addr.330, align 8, !tbaa !1403, !dbg !1603
	br label %L..inline.10329
L..inline.10329:
	%109 = load i32, i32* %.ndi0006.addr, align 4, !tbaa !1401, !dbg !1603
	%110 = sext i32  %109 to i64, !dbg !1603
	%111 = load double*, double** %..inline.addr.8, align 8, !tbaa !1398, !dbg !1603
	%112 = getelementptr double, double*  %111, i64  %110, !dbg !1603
	%113 = load double, double*  %112, align 8, !tbaa !1399, !dbg !1603
	%114 = bitcast [20 x double]* %..inline.addr.370 to double*, !dbg !1603
	%115 = getelementptr double, double*  %114, i64  %110, !dbg !1603
	store double  %113, double*  %115, align 8, !tbaa !1399, !dbg !1603
	%116 = bitcast [20 x double]* %..inline.addr.380 to double*, !dbg !1603
	%117 = getelementptr double, double*  %116, i64  %110, !dbg !1603
	store double  %113, double*  %117, align 8, !tbaa !1399, !dbg !1603
	%118 = add i32  %109, 1, !dbg !1603
	store i32  %118, i32* %.ndi0006.addr, align 4, !tbaa !1401, !dbg !1603
	%119 = load i32, i32* %.lcr013074.addr, align 4, !tbaa !1401, !dbg !1603
	%120 = icmp slt i32  %118,  %119, !dbg !1603
	br i1  %120, label %L.B0197, label %L.B0198, !dbg !1603
L.B0198:
	%121 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1605
	%122 = sext i32  %121 to i64, !dbg !1605
	%123 = load double*, double** %..inline.addr.8, align 8, !tbaa !1398, !dbg !1605
	%124 = getelementptr double, double*  %123, i64  %122, !dbg !1605
	%125 = load double, double*  %124, align 8, !tbaa !1399, !dbg !1605
	store double  %125, double* %result.addr, align 8, !tbaa !1403, !dbg !1605

	%126 = sub i32  %121, 1, !dbg !1607
	store i32  %126, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1607
	%127 = load i32, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !1598
	%128 = sext i32  %127 to i64, !dbg !1598
	%129 = sub i64  %128, 1, !dbg !1598
	store i64  %129, i64* %.TRP0002.addr, align 8, !tbaa !1449, !dbg !1598
	store i64 0, i64* %.ndk0012.addr, align 8, !tbaa !1451, !dbg !1603
	%130 = icmp sle i64  %129, 0, !dbg !1603
	br i1  %130, label %L.B0200, label %L.B0208, !dbg !1603
L.B0208:
	store i64  %129, i64* %.lcr053074.addr, align 8, !tbaa !1451, !dbg !1603
	br label %L.B0199
L.B0199:
	%131 = load i32, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !1603
	%132 = load i64, i64* %.ndk0012.addr, align 8, !tbaa !1451, !dbg !1603
	%133 = trunc i64  %132 to i32, !dbg !1603
	%134 = sub i32  %131,  %133, !dbg !1603
	%135 = sext i32  %134 to i64, !dbg !1603
	%136 = sub i64  %135, 1, !dbg !1603
	store i64  %136, i64* %.TRP0003.addr, align 8, !tbaa !1449, !dbg !1603
	store i64 0, i64* %.ndk0014.addr, align 8, !tbaa !1451, !dbg !1603
	%137 = icmp sle i64  %136, 0, !dbg !1603
	br i1  %137, label %L.B0202, label %L.B0209, !dbg !1603
L.B0209:
	store i64  %136, i64* %.lcr053075.addr, align 8, !tbaa !1451, !dbg !1608
	store i32  %133, i32* %.r4.0386.addr, align 4, !tbaa !1401, !dbg !1608
	%138 = bitcast [20 x double]* %..inline.addr.380 to i8*, !dbg !1608
	store i8*  %138, i8** %.G0020.addr, align 8, !tbaa !1398, !dbg !1608
	br label %L.B0201
L.B0201:
	%139 = load i64, i64* %.ndk0014.addr, align 8, !tbaa !1451, !dbg !1603
	%140 = trunc i64  %139 to i32, !dbg !1603
	%141 = sext i32  %140 to i64, !dbg !1603
	%142 = load double*, double** %..inline.addr.7, align 8, !tbaa !1398, !dbg !1603
	%143 = getelementptr double, double*  %142, i64  %141, !dbg !1603
	%144 = load double, double*  %143, align 8, !tbaa !1399, !dbg !1603
	%145 = load i8*, i8** %.G0020.addr, align 8, !tbaa !1398, !dbg !1603
	%146 = getelementptr i8, i8*  %145, i64 8, !dbg !1603
	%147 = bitcast i8*  %146 to double*, !dbg !1603
	%148 = load double, double*  %147, align 8, !tbaa !1399, !dbg !1603
	%149 = bitcast [20 x double]* %..inline.addr.370 to double*, !dbg !1603
	%150 = getelementptr double, double*  %149, i64  %141, !dbg !1603
	%151 = load double, double*  %150, align 8, !tbaa !1399, !dbg !1603
	%152 = fsub double  %148,  %151, !dbg !1603
	%153 = load i32, i32* %.r4.0386.addr, align 4, !tbaa !1401, !dbg !1603
	%154 = add i32  %153,  %140, !dbg !1603
	%155 = sext i32  %154 to i64, !dbg !1603
	%156 = bitcast double*  %142 to i8*, !dbg !1603
	%157 = getelementptr i8, i8*  %156, i64 8, !dbg !1603
	%158 = bitcast i8*  %157 to double*, !dbg !1603
	%159 = getelementptr double, double*  %158, i64  %155, !dbg !1603
	%160 = load double, double*  %159, align 8, !tbaa !1399, !dbg !1603
	%161 = fsub double  %144,  %160, !dbg !1603
	%162 = fdiv double  %152,  %161, !dbg !1603
	%163 = fmul double  %162,  %160, !dbg !1603
	store double  %163, double*  %150, align 8, !tbaa !1399, !dbg !1603
	%164 = fmul double  %144,  %162, !dbg !1603
	%165 = bitcast i8*  %145 to double*, !dbg !1603
	store double  %164, double*  %165, align 8, !tbaa !1399, !dbg !1603
	%166 = add i64  %139, 1, !dbg !1603
	store i64  %166, i64* %.ndk0014.addr, align 8, !tbaa !1451, !dbg !1603
	store i8*  %146, i8** %.G0020.addr, align 8, !tbaa !1398, !dbg !1608
	%167 = load i64, i64* %.lcr053075.addr, align 8, !tbaa !1451, !dbg !1603
	%168 = icmp slt i64  %166,  %167, !dbg !1603
	br i1  %168, label %L.B0201, label %L.B0202, !dbg !1603
L.B0202:
	%169 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1603
	%170 = add i32  %169, 1, !dbg !1603
	%171 = mul i32  %170, 2, !dbg !1603
	%172 = load i32, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !1603
	%173 = load i64, i64* %.ndk0012.addr, align 8, !tbaa !1451, !dbg !1603
	%174 = trunc i64  %173 to i32, !dbg !1603
	%175 = sub i32  %172,  %174, !dbg !1603
	%176 = sub i32  %175, 1, !dbg !1603
	%177 = icmp sge i32  %171,  %176, !dbg !1603
	br i1  %177, label %L..inline.10349, label %L.B0210, !dbg !1603
L.B0210:
	%178 = sext i32  %169 to i64, !dbg !1603
	%179 = bitcast [20 x double]* %..inline.addr.380 to i8*, !dbg !1603
	%180 = getelementptr i8, i8*  %179, i64 8, !dbg !1603
	%181 = bitcast i8*  %180 to double*, !dbg !1603
	%182 = getelementptr double, double*  %181, i64  %178, !dbg !1603
	%183 = load double, double*  %182, align 8, !tbaa !1399, !dbg !1603
	store double  %183, double* %dresult.addr, align 8, !tbaa !1403, !dbg !1603
	br label %L..inline.10351, !dbg !1603
L..inline.10349:
	%184 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1603
	%185 = sext i32  %184 to i64, !dbg !1603
	%186 = bitcast [20 x double]* %..inline.addr.370 to double*, !dbg !1603
	%187 = getelementptr double, double*  %186, i64  %185, !dbg !1603
	%188 = load double, double*  %187, align 8, !tbaa !1399, !dbg !1603
	store double  %188, double* %dresult.addr, align 8, !tbaa !1403, !dbg !1603

	%189 = sub i32  %184, 1, !dbg !1603
	store i32  %189, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1603
	br label %L..inline.10351
L..inline.10351:
	%190 = load double, double* %dresult.addr, align 8, !tbaa !1403, !dbg !1603
	%191 = load double, double* %result.addr, align 8, !tbaa !1403, !dbg !1603
	%192 = fadd double  %190,  %191, !dbg !1603
	store double  %192, double* %result.addr, align 8, !tbaa !1403, !dbg !1603
	%193 = load i64, i64* %.ndk0012.addr, align 8, !tbaa !1451, !dbg !1603
	%194 = add i64  %193, 1, !dbg !1603
	store i64  %194, i64* %.ndk0012.addr, align 8, !tbaa !1451, !dbg !1603
	%195 = load i64, i64* %.lcr053074.addr, align 8, !tbaa !1451, !dbg !1603
	%196 = icmp slt i64  %194,  %195, !dbg !1603
	br i1  %196, label %L.B0199, label %L.B0200, !dbg !1603
L.B0200:
	%197 = load double, double* %dresult.addr, align 8, !tbaa !1403, !dbg !1609
	%198 = call double @llvm.fabs.f64 (double  %197), !dbg !1609
	%199 = load double, double* %absacc.addr, align 8, !tbaa !1403, !dbg !1609
	%200 = fcmp olt double  %198,  %199, !dbg !1609
	br i1  %200, label %L.N0005, label %L.B0035, !dbg !1609
L.B0035:
	%201 = load i8*, i8** %.G0035.addr, align 8, !tbaa !1398, !dbg !1610
	%202 = getelementptr i8, i8*  %201, i64 18446744073709551608, !dbg !1610
	%203 = bitcast i8*  %202 to double*, !dbg !1610
	%204 = load double, double*  %203, align 8, !tbaa !1399, !dbg !1610
	%205 = bitcast i8*  %201 to double*, !dbg !1610
	store double  %204, double*  %205, align 8, !tbaa !1399, !dbg !1610
	%206 = load i8*, i8** %.G0033.addr, align 8, !tbaa !1398, !dbg !1611
	%207 = getelementptr i8, i8*  %206, i64 18446744073709551608, !dbg !1611
	%208 = bitcast i8*  %207 to double*, !dbg !1611
	%209 = load double, double*  %208, align 8, !tbaa !1399, !dbg !1611
	%210 = fmul double  %209,  2.50000000000000000E-1, !dbg !1611
	%211 = bitcast i8*  %206 to double*, !dbg !1611
	store double  %210, double*  %211, align 8, !tbaa !1399, !dbg !1611
	%212 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1612

	%213 = add i32  %212, 1, !dbg !1612
	store i32  %213, i32* %j.addr, align 4, !tbaa !1401, !dbg !1612
	%214 = getelementptr i8, i8*  %201, i64 8, !dbg !1582
	store i8*  %214, i8** %.G0035.addr, align 8, !tbaa !1398, !dbg !1582
	%215 = getelementptr i8, i8*  %206, i64 8, !dbg !1582
	store i8*  %215, i8** %.G0033.addr, align 8, !tbaa !1398, !dbg !1582
	%216 = load i32, i32* %.x3073.addr, align 4, !tbaa !1401, !dbg !1582
	%217 = sub i32  %216, 1, !dbg !1582
	store i32  %217, i32* %.x3073.addr, align 4, !tbaa !1401, !dbg !1582
	%218 = icmp sgt i32  %217, 0, !dbg !1612
	br i1  %218, label %L.B0033, label %L.N0007, !dbg !1612
L.N0007:
	br label %L.N0005
L.N0005:
	%219 = load double, double* %result.addr, align 8, !tbaa !1403, !dbg !1613
	ret double  %219, !dbg !1613
L.N0006:
	ret double 0.0
}
define double @_Z7RombergRK11T1DFunctionddd(%struct.T1DFunction* %func.arg, double %a.arg, double %b.arg, double %absacc.arg) #0 !dbg !1616 {
L.entry:
	%func.addr = alloca %struct.T1DFunction*, align 8
	%a.addr = alloca double, align 8
	%b.addr = alloca double, align 8
	%absacc.addr = alloca double, align 8
	%it.addr = alloca i32, align 4
	%H.addr = alloca [9 x double], align 8
	%dresult_pol.addr = alloca double, align 8
	%j.addr = alloca i32, align 4
	%.x4097.addr = alloca i32, align 4
	%S.addr = alloca [9 x double], align 8
	%.G0052.addr = alloca i8*, align 8
	%.G0050.addr = alloca i8*, align 8
	%..inline.addr = alloca %struct.T1DFunction*, align 8
	%..inline.addr.1 = alloca double, align 8
	%..inline.addr.2 = alloca double, align 8
	%..inline.addr.3 = alloca double*, align 8
	%.Q0007.addr = alloca double, align 8
	%.Q0008.addr = alloca double, align 8
	%..inline.addr.4 = alloca double, align 8
	%..inline.addr.5 = alloca double, align 8
	%..inline.addr.6 = alloca double, align 8
	%.x4098.addr = alloca i32, align 4
	%.Q0009.addr = alloca double, align 8
	%result.addr = alloca double, align 8
	%k1.addr = alloca i32, align 4
	%..inline.addr.7 = alloca double*, align 8
	%..inline.addr.8 = alloca double*, align 8
	%..inline.addr.9 = alloca i32, align 4
	%..inline.addr.320 = alloca i32, align 4
	%..inline.addr.330 = alloca double, align 8
	%.ndi0010.addr = alloca i32, align 4
	%.lcr014098.addr = alloca i32, align 4
	%..inline.addr.360 = alloca double, align 8
	%..inline.addr.370 = alloca [20 x double], align 8
	%..inline.addr.380 = alloca [20 x double], align 8
	%result_pol.addr = alloca double, align 8
	%.TRP0004.addr = alloca i64, align 8
	%.ndk0020.addr = alloca i64, align 8
	%.lcr054098.addr = alloca i64, align 8
	%.TRP0005.addr = alloca i64, align 8
	%.ndk0022.addr = alloca i64, align 8
	%.lcr054099.addr = alloca i64, align 8
	%.r5.0631.addr = alloca i32, align 4
	%.G0037.addr = alloca i8*, align 8
	%..inline.addr.480 = alloca double, align 8
	%..inline.addr.490 = alloca double, align 8
	%..inline.addr.500 = alloca double, align 8
	%..inline.addr.510 = alloca double*, align 8
	%..inline.addr.520 = alloca double*, align 8
	%..inline.addr.530 = alloca i32, align 4
	%..inline.addr.540 = alloca i32, align 4
	%..inline.addr.550 = alloca double, align 8
	%..inline.addr.560 = alloca i32, align 4
	%..inline.addr.570 = alloca double, align 8
	%result_rat.addr = alloca double, align 8
	%dresult_rat.addr = alloca double, align 8
	%..inline.addr.600 = alloca [20 x double], align 8
	%..inline.addr.610 = alloca [20 x double], align 8
	%.D0019.addr = alloca i32, align 4
	%..inline.addr.630 = alloca i32, align 4
	%.lcr014099.addr = alloca i32, align 4
	%..inline.addr.650 = alloca double, align 8
	%..inline.addr.660 = alloca double, align 8
	%..inline.addr.670 = alloca double, align 8
	%..inline.addr.680 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.690 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.700 = alloca i32, align 4
	%..inline.addr.710 = alloca i64, align 8
	%.Q0010.addr = alloca %struct._ZSo*, align 8
	%..inline.addr.730 = alloca double, align 8
	%.D0020.addr = alloca i32, align 4
	%.I0003.addr = alloca double, align 8
	%dresult_polrat.addr = alloca double, align 8
	%dresult.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %func.addr, metadata !1668, metadata !1366), !dbg !1617
	store %struct.T1DFunction* %func.arg, %struct.T1DFunction** %func.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %func.addr, metadata !1669, metadata !1366), !dbg !1617
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1670, metadata !1366), !dbg !1617
	store double %a.arg, double* %a.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1671, metadata !1366), !dbg !1617
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1672, metadata !1366), !dbg !1617
	store double %b.arg, double* %b.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1673, metadata !1366), !dbg !1617
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !1674, metadata !1366), !dbg !1617
	store double %absacc.arg, double* %absacc.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !1675, metadata !1366), !dbg !1617
	call void @llvm.dbg.declare (metadata i32* %it.addr, metadata !1677, metadata !1366), !dbg !1617
	store i32 0, i32* %it.addr, align 4, !tbaa !1401, !dbg !1676
	call void @llvm.dbg.declare (metadata [9 x double]* %H.addr, metadata !1679, metadata !1366), !dbg !1617
	%0 = bitcast [9 x double]* %H.addr to double*, !dbg !1678
	store double  1.00000000000000000E+0, double*  %0, align 8, !tbaa !1399, !dbg !1678
	call void @llvm.dbg.declare (metadata double* %dresult_pol.addr, metadata !1681, metadata !1366), !dbg !1617
	store double  0.00000000000000000E+0, double* %dresult_pol.addr, align 8, !tbaa !1403, !dbg !1680
	call void @llvm.dbg.declare (metadata i32* %j.addr, metadata !1683, metadata !1366), !dbg !1617
	store i32 1, i32* %j.addr, align 4, !tbaa !1401, !dbg !1682
	store i32 8, i32* %.x4097.addr, align 4, !tbaa !1401, !dbg !1682
	call void @llvm.dbg.declare (metadata [9 x double]* %S.addr, metadata !1684, metadata !1366), !dbg !1617
	%1 = bitcast [9 x double]* %S.addr to i8*, !dbg !1682
	%2 = getelementptr i8, i8*  %1, i64 8, !dbg !1682
	store i8*  %2, i8** %.G0052.addr, align 8, !tbaa !1398, !dbg !1682
	%3 = bitcast [9 x double]* %H.addr to i8*, !dbg !1682
	%4 = getelementptr i8, i8*  %3, i64 8, !dbg !1682
	store i8*  %4, i8** %.G0050.addr, align 8, !tbaa !1398, !dbg !1682
	br label %L.B0039
L.B0039:
	%5 = load %struct.T1DFunction*, %struct.T1DFunction** %func.addr, align 8, !tbaa !1398, !dbg !1685
	%6 = bitcast %struct.T1DFunction*  %5 to i8*, !dbg !1685
	%7 = bitcast %struct.T1DFunction** %..inline.addr to i8**, !dbg !1685
	store i8*  %6, i8**  %7, align 8, !tbaa !1398, !dbg !1685
	%8 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !1685
	store double  %8, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1685
	%9 = load double, double* %b.addr, align 8, !tbaa !1403, !dbg !1685
	store double  %9, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !1685
	%10 = bitcast [9 x double]* %S.addr to i8*, !dbg !1685
	%11 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1685
	%12 = sext i32  %11 to i64, !dbg !1685
	%13 = sub i64  %12, 1, !dbg !1685
	%14 = mul i64  %13, 8, !dbg !1685
	%15 = getelementptr i8, i8*  %10, i64  %14, !dbg !1685
	%16 = bitcast double** %..inline.addr.3 to i8**, !dbg !1685
	store i8*  %15, i8**  %16, align 8, !tbaa !1398, !dbg !1685
	%17 = icmp ne i32  %11, 1, !dbg !1685
	br i1  %17, label %L..inline.10790, label %L.B0232, !dbg !1685
L.B0232:
	%18 = load %struct.T1DFunction*, %struct.T1DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !1686
	%19 = bitcast %struct.T1DFunction*  %18 to double (%struct.T1DFunction*, double)***, !dbg !1686
	%20 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %19, align 8, !tbaa !1398, !dbg !1686
	%21 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %20, align 8, !tbaa !1398, !dbg !1686
	%22 = call double  %21 (%struct.T1DFunction*  %18, double  %9), !dbg !1686
	store double  %22, double* %.Q0007.addr, align 8, !tbaa !1403, !dbg !1686
	%23 = load %struct.T1DFunction*, %struct.T1DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !1686
	%24 = load double, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1686
	%25 = bitcast %struct.T1DFunction*  %23 to double (%struct.T1DFunction*, double)***, !dbg !1686
	%26 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %25, align 8, !tbaa !1398, !dbg !1686
	%27 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %26, align 8, !tbaa !1398, !dbg !1686
	%28 = call double  %27 (%struct.T1DFunction*  %23, double  %24), !dbg !1686
	store double  %28, double* %.Q0008.addr, align 8, !tbaa !1403, !dbg !1686
	%29 = load double, double* %.Q0007.addr, align 8, !tbaa !1403, !dbg !1686
	%30 = fadd double  %29,  %28, !dbg !1686
	%31 = load double, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !1686
	%32 = load double, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1686
	%33 = fsub double  %31,  %32, !dbg !1686
	%34 = fmul double  %33,  5.00000000000000000E-1, !dbg !1686
	%35 = fmul double  %30,  %34, !dbg !1686
	%36 = load double*, double** %..inline.addr.3, align 8, !tbaa !1398, !dbg !1686
	store double  %35, double*  %36, align 8, !tbaa !1399, !dbg !1686
	store i32 1, i32* %it.addr, align 4, !tbaa !1401, !dbg !1686
	br label %L..inline.10796, !dbg !1687
L..inline.10790:
	%37 = load double, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !1688
	%38 = load double, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1688
	%39 = fsub double  %37,  %38, !dbg !1688
	%40 = load i32, i32* %it.addr, align 4, !tbaa !1401, !dbg !1688
	%41 = sitofp i32  %40 to double, !dbg !1688
	%42 = fdiv double  %39,  %41, !dbg !1688
	store double  %42, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !1688
	%43 = call double @llvm.fma.f64 (double  %42, double  5.00000000000000000E-1, double  %38), !dbg !1689
	store double  %43, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1689
	store double  0.00000000000000000E+0, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1690
	%44 = icmp sle i32  %40, 0, !dbg !1691
	br i1  %44, label %L..inline.10805, label %L.B0233, !dbg !1691
L.B0233:
	store i32  %40, i32* %.x4098.addr, align 4, !tbaa !1401, !dbg !1691
	br label %L..inline.10804
L..inline.10804:
	%45 = load %struct.T1DFunction*, %struct.T1DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !1692
	%46 = load double, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1692
	%47 = bitcast %struct.T1DFunction*  %45 to double (%struct.T1DFunction*, double)***, !dbg !1692
	%48 = load double (%struct.T1DFunction*, double)**, double (%struct.T1DFunction*, double)***  %47, align 8, !tbaa !1398, !dbg !1692
	%49 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %48, align 8, !tbaa !1398, !dbg !1692
	%50 = call double  %49 (%struct.T1DFunction*  %45, double  %46), !dbg !1692
	store double  %50, double* %.Q0009.addr, align 8, !tbaa !1403, !dbg !1692
	%51 = load double, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1692
	%52 = fadd double  %50,  %51, !dbg !1692
	store double  %52, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1692
	%53 = load double, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !1692
	%54 = load double, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1692
	%55 = fadd double  %53,  %54, !dbg !1692
	store double  %55, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1692
	%56 = load i32, i32* %.x4098.addr, align 4, !tbaa !1401, !dbg !1691
	%57 = sub i32  %56, 1, !dbg !1691
	store i32  %57, i32* %.x4098.addr, align 4, !tbaa !1401, !dbg !1691
	%58 = icmp sgt i32  %57, 0, !dbg !1693
	br i1  %58, label %L..inline.10804, label %L..inline.10805, !dbg !1693
L..inline.10805:
	%59 = load double, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !1693
	%60 = load double, double* %..inline.addr.1, align 8, !tbaa !1403, !dbg !1693
	%61 = fsub double  %59,  %60, !dbg !1693
	%62 = load double, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1693
	%63 = fmul double  %61,  %62, !dbg !1693
	%64 = load i32, i32* %it.addr, align 4, !tbaa !1401, !dbg !1693
	%65 = sitofp i32  %64 to double, !dbg !1693
	%66 = fdiv double  %63,  %65, !dbg !1693
	%67 = load double*, double** %..inline.addr.3, align 8, !tbaa !1398, !dbg !1693
	%68 = load double, double*  %67, align 8, !tbaa !1399, !dbg !1693
	%69 = fadd double  %66,  %68, !dbg !1693
	%70 = fmul double  %69,  5.00000000000000000E-1, !dbg !1693
	store double  %70, double*  %67, align 8, !tbaa !1399, !dbg !1693
	%71 = mul i32  %64, 2, !dbg !1694
	store i32  %71, i32* %it.addr, align 4, !tbaa !1401, !dbg !1694
	br label %L..inline.10796
L..inline.10796:
	%72 = load i8*, i8** %.G0052.addr, align 8, !tbaa !1398, !dbg !1695
	%73 = getelementptr i8, i8*  %72, i64 18446744073709551608, !dbg !1695
	%74 = bitcast i8*  %73 to double*, !dbg !1695
	%75 = load double, double*  %74, align 8, !tbaa !1399, !dbg !1695
	call void @llvm.dbg.declare (metadata double* %result.addr, metadata !1696, metadata !1366), !dbg !1617
	store double  %75, double* %result.addr, align 8, !tbaa !1403, !dbg !1695
	%76 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1697
	%77 = icmp slt i32  %76, 2, !dbg !1697
	br i1  %77, label %L.B0041, label %L.B0234, !dbg !1697
L.B0234:
	%78 = icmp slt i32  %76, 5, !dbg !1698
	%79 = select i1  %78, i32  %76, i32 5, !dbg !1698
	call void @llvm.dbg.declare (metadata i32* %k1.addr, metadata !1699, metadata !1366), !dbg !1617
	store i32  %79, i32* %k1.addr, align 4, !tbaa !1401, !dbg !1698
	%80 = bitcast [9 x double]* %H.addr to i8*, !dbg !1700
	%81 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1700
	%82 = sub i32  %81,  %79, !dbg !1700
	%83 = sext i32  %82 to i64, !dbg !1700
	%84 = mul i64  %83, 8, !dbg !1700
	%85 = getelementptr i8, i8*  %80, i64  %84, !dbg !1700
	%86 = bitcast double** %..inline.addr.7 to i8**, !dbg !1700
	store i8*  %85, i8**  %86, align 8, !tbaa !1398, !dbg !1700
	%87 = bitcast [9 x double]* %S.addr to i8*, !dbg !1700
	%88 = mul i64  %83, 8, !dbg !1700
	%89 = getelementptr i8, i8*  %87, i64  %88, !dbg !1700
	%90 = bitcast double** %..inline.addr.8 to i8**, !dbg !1700
	store i8*  %89, i8**  %90, align 8, !tbaa !1398, !dbg !1700
	store i32  %79, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !1700
	store i32 0, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1701
	%91 = load double*, double** %..inline.addr.7, align 8, !tbaa !1398, !dbg !1702
	%92 = load double, double*  %91, align 8, !tbaa !1399, !dbg !1702
	%93 = fsub double -0.00000000e+00,  %92, !dbg !1702
	%94 = call double @llvm.fabs.f64 (double  %93), !dbg !1702
	store double  %94, double* %..inline.addr.330, align 8, !tbaa !1403, !dbg !1702
	store i32 0, i32* %.ndi0010.addr, align 4, !tbaa !1401, !dbg !1703
	%95 = icmp slt i32  %81, 5, !dbg !1703
	%96 = select i1  %95, i32  %81, i32 5, !dbg !1703
	%97 = icmp sle i32  %96, 0, !dbg !1703
	br i1  %97, label %L.B0227, label %L.B0235, !dbg !1703
L.B0235:
	%98 = icmp slt i32  %81, 5, !dbg !1704
	%99 = select i1  %98, i32  %81, i32 5, !dbg !1704
	store i32  %99, i32* %.lcr014098.addr, align 4, !tbaa !1401, !dbg !1704
	br label %L.B0226
L.B0226:
	%100 = load i32, i32* %.ndi0010.addr, align 4, !tbaa !1401, !dbg !1703
	%101 = sext i32  %100 to i64, !dbg !1703
	%102 = load double*, double** %..inline.addr.7, align 8, !tbaa !1398, !dbg !1703
	%103 = getelementptr double, double*  %102, i64  %101, !dbg !1703
	%104 = load double, double*  %103, align 8, !tbaa !1399, !dbg !1703
	%105 = fsub double -0.00000000e+00,  %104, !dbg !1703
	%106 = call double @llvm.fabs.f64 (double  %105), !dbg !1703
	store double  %106, double* %..inline.addr.360, align 8, !tbaa !1403, !dbg !1703
	%107 = load double, double* %..inline.addr.330, align 8, !tbaa !1403, !dbg !1703
	%108 = fcmp uge double  %106,  %107, !dbg !1703
	br i1  %108, label %L..inline.10823, label %L.B0236, !dbg !1703
L.B0236:
	store i32  %100, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1703
	store double  %106, double* %..inline.addr.330, align 8, !tbaa !1403, !dbg !1703
	br label %L..inline.10823
L..inline.10823:
	%109 = load i32, i32* %.ndi0010.addr, align 4, !tbaa !1401, !dbg !1703
	%110 = sext i32  %109 to i64, !dbg !1703
	%111 = load double*, double** %..inline.addr.8, align 8, !tbaa !1398, !dbg !1703
	%112 = getelementptr double, double*  %111, i64  %110, !dbg !1703
	%113 = load double, double*  %112, align 8, !tbaa !1399, !dbg !1703
	%114 = bitcast [20 x double]* %..inline.addr.370 to double*, !dbg !1703
	%115 = getelementptr double, double*  %114, i64  %110, !dbg !1703
	store double  %113, double*  %115, align 8, !tbaa !1399, !dbg !1703
	%116 = bitcast [20 x double]* %..inline.addr.380 to double*, !dbg !1703
	%117 = getelementptr double, double*  %116, i64  %110, !dbg !1703
	store double  %113, double*  %117, align 8, !tbaa !1399, !dbg !1703
	%118 = add i32  %109, 1, !dbg !1703
	store i32  %118, i32* %.ndi0010.addr, align 4, !tbaa !1401, !dbg !1703
	%119 = load i32, i32* %.lcr014098.addr, align 4, !tbaa !1401, !dbg !1703
	%120 = icmp slt i32  %118,  %119, !dbg !1703
	br i1  %120, label %L.B0226, label %L.B0227, !dbg !1703
L.B0227:
	%121 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1705
	%122 = sext i32  %121 to i64, !dbg !1705
	%123 = load double*, double** %..inline.addr.8, align 8, !tbaa !1398, !dbg !1705
	%124 = getelementptr double, double*  %123, i64  %122, !dbg !1705
	%125 = load double, double*  %124, align 8, !tbaa !1399, !dbg !1705
	call void @llvm.dbg.declare (metadata double* %result_pol.addr, metadata !1706, metadata !1366), !dbg !1617
	store double  %125, double* %result_pol.addr, align 8, !tbaa !1403, !dbg !1705
	%126 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1707

	%127 = sub i32  %126, 1, !dbg !1708
	store i32  %127, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1708
	%128 = load i32, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !1709
	%129 = sext i32  %128 to i64, !dbg !1709
	%130 = sub i64  %129, 1, !dbg !1709
	store i64  %130, i64* %.TRP0004.addr, align 8, !tbaa !1449, !dbg !1709
	store i64 0, i64* %.ndk0020.addr, align 8, !tbaa !1451, !dbg !1703
	%131 = icmp sle i64  %130, 0, !dbg !1703
	br i1  %131, label %L.B0229, label %L.B0237, !dbg !1703
L.B0237:
	store i64  %130, i64* %.lcr054098.addr, align 8, !tbaa !1451, !dbg !1710
	br label %L.B0228
L.B0228:
	%132 = load i32, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !1703
	%133 = load i64, i64* %.ndk0020.addr, align 8, !tbaa !1451, !dbg !1703
	%134 = trunc i64  %133 to i32, !dbg !1703
	%135 = sub i32  %132,  %134, !dbg !1703
	%136 = sext i32  %135 to i64, !dbg !1703
	%137 = sub i64  %136, 1, !dbg !1703
	store i64  %137, i64* %.TRP0005.addr, align 8, !tbaa !1449, !dbg !1703
	store i64 0, i64* %.ndk0022.addr, align 8, !tbaa !1451, !dbg !1703
	%138 = icmp sle i64  %137, 0, !dbg !1703
	br i1  %138, label %L.B0231, label %L.B0238, !dbg !1703
L.B0238:
	store i64  %137, i64* %.lcr054099.addr, align 8, !tbaa !1451, !dbg !1711
	store i32  %134, i32* %.r5.0631.addr, align 4, !tbaa !1401, !dbg !1711
	%139 = bitcast [20 x double]* %..inline.addr.380 to i8*, !dbg !1711
	store i8*  %139, i8** %.G0037.addr, align 8, !tbaa !1398, !dbg !1711
	br label %L.B0230
L.B0230:
	%140 = load i64, i64* %.ndk0022.addr, align 8, !tbaa !1451, !dbg !1703
	%141 = trunc i64  %140 to i32, !dbg !1703
	%142 = sext i32  %141 to i64, !dbg !1703
	%143 = load double*, double** %..inline.addr.7, align 8, !tbaa !1398, !dbg !1703
	%144 = getelementptr double, double*  %143, i64  %142, !dbg !1703
	%145 = load double, double*  %144, align 8, !tbaa !1399, !dbg !1703
	%146 = load i8*, i8** %.G0037.addr, align 8, !tbaa !1398, !dbg !1703
	%147 = getelementptr i8, i8*  %146, i64 8, !dbg !1703
	%148 = bitcast i8*  %147 to double*, !dbg !1703
	%149 = load double, double*  %148, align 8, !tbaa !1399, !dbg !1703
	%150 = bitcast [20 x double]* %..inline.addr.370 to double*, !dbg !1703
	%151 = getelementptr double, double*  %150, i64  %142, !dbg !1703
	%152 = load double, double*  %151, align 8, !tbaa !1399, !dbg !1703
	%153 = fsub double  %149,  %152, !dbg !1703
	%154 = load i32, i32* %.r5.0631.addr, align 4, !tbaa !1401, !dbg !1703
	%155 = add i32  %154,  %141, !dbg !1703
	%156 = sext i32  %155 to i64, !dbg !1703
	%157 = bitcast double*  %143 to i8*, !dbg !1703
	%158 = getelementptr i8, i8*  %157, i64 8, !dbg !1703
	%159 = bitcast i8*  %158 to double*, !dbg !1703
	%160 = getelementptr double, double*  %159, i64  %156, !dbg !1703
	%161 = load double, double*  %160, align 8, !tbaa !1399, !dbg !1703
	%162 = fsub double  %145,  %161, !dbg !1703
	%163 = fdiv double  %153,  %162, !dbg !1703
	%164 = fmul double  %163,  %161, !dbg !1703
	store double  %164, double*  %151, align 8, !tbaa !1399, !dbg !1703
	%165 = fmul double  %145,  %163, !dbg !1703
	%166 = bitcast i8*  %146 to double*, !dbg !1703
	store double  %165, double*  %166, align 8, !tbaa !1399, !dbg !1703
	%167 = add i64  %140, 1, !dbg !1703
	store i64  %167, i64* %.ndk0022.addr, align 8, !tbaa !1451, !dbg !1703
	store i8*  %147, i8** %.G0037.addr, align 8, !tbaa !1398, !dbg !1711
	%168 = load i64, i64* %.lcr054099.addr, align 8, !tbaa !1451, !dbg !1703
	%169 = icmp slt i64  %167,  %168, !dbg !1703
	br i1  %169, label %L.B0230, label %L.B0231, !dbg !1703
L.B0231:
	%170 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1703
	%171 = add i32  %170, 1, !dbg !1703
	%172 = mul i32  %171, 2, !dbg !1703
	%173 = load i32, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !1703
	%174 = load i64, i64* %.ndk0020.addr, align 8, !tbaa !1451, !dbg !1703
	%175 = trunc i64  %174 to i32, !dbg !1703
	%176 = sub i32  %173,  %175, !dbg !1703
	%177 = sub i32  %176, 1, !dbg !1703
	%178 = icmp sge i32  %172,  %177, !dbg !1703
	br i1  %178, label %L..inline.10843, label %L.B0239, !dbg !1703
L.B0239:
	%179 = sext i32  %170 to i64, !dbg !1703
	%180 = bitcast [20 x double]* %..inline.addr.380 to i8*, !dbg !1703
	%181 = getelementptr i8, i8*  %180, i64 8, !dbg !1703
	%182 = bitcast i8*  %181 to double*, !dbg !1703
	%183 = getelementptr double, double*  %182, i64  %179, !dbg !1703
	%184 = load double, double*  %183, align 8, !tbaa !1399, !dbg !1703
	store double  %184, double* %dresult_pol.addr, align 8, !tbaa !1403, !dbg !1703
	br label %L..inline.10845, !dbg !1703
L..inline.10843:
	%185 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1703
	%186 = sext i32  %185 to i64, !dbg !1703
	%187 = bitcast [20 x double]* %..inline.addr.370 to double*, !dbg !1703
	%188 = getelementptr double, double*  %187, i64  %186, !dbg !1703
	%189 = load double, double*  %188, align 8, !tbaa !1399, !dbg !1703
	store double  %189, double* %dresult_pol.addr, align 8, !tbaa !1403, !dbg !1703

	%190 = sub i32  %185, 1, !dbg !1703
	store i32  %190, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !1703
	br label %L..inline.10845
L..inline.10845:
	%191 = load double, double* %dresult_pol.addr, align 8, !tbaa !1403, !dbg !1703
	%192 = load double, double* %result_pol.addr, align 8, !tbaa !1403, !dbg !1703
	%193 = fadd double  %191,  %192, !dbg !1703
	store double  %193, double* %result_pol.addr, align 8, !tbaa !1403, !dbg !1703
	%194 = load i64, i64* %.ndk0020.addr, align 8, !tbaa !1451, !dbg !1703
	%195 = add i64  %194, 1, !dbg !1703
	store i64  %195, i64* %.ndk0020.addr, align 8, !tbaa !1451, !dbg !1703
	%196 = load i64, i64* %.lcr054098.addr, align 8, !tbaa !1451, !dbg !1703
	%197 = icmp slt i64  %195,  %196, !dbg !1703
	br i1  %197, label %L.B0228, label %L.B0229, !dbg !1703
L.B0229:
	%198 = bitcast [9 x double]* %H.addr to i8*, !dbg !1712
	%199 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1712
	%200 = load i32, i32* %k1.addr, align 4, !tbaa !1401, !dbg !1712
	%201 = sub i32  %199,  %200, !dbg !1712
	%202 = sext i32  %201 to i64, !dbg !1712
	%203 = mul i64  %202, 8, !dbg !1712
	%204 = getelementptr i8, i8*  %198, i64  %203, !dbg !1712
	%205 = bitcast double** %..inline.addr.510 to i8**, !dbg !1712
	store i8*  %204, i8**  %205, align 8, !tbaa !1398, !dbg !1712
	%206 = bitcast [9 x double]* %S.addr to i8*, !dbg !1712
	%207 = mul i64  %202, 8, !dbg !1712
	%208 = getelementptr i8, i8*  %206, i64  %207, !dbg !1712
	%209 = bitcast double** %..inline.addr.520 to i8**, !dbg !1712
	store i8*  %208, i8**  %209, align 8, !tbaa !1398, !dbg !1712
	store i32  %200, i32* %..inline.addr.530, align 4, !tbaa !1401, !dbg !1712
	store i32 0, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !1712
	%210 = load double*, double** %..inline.addr.510, align 8, !tbaa !1398, !dbg !1713
	%211 = load double, double*  %210, align 8, !tbaa !1399, !dbg !1713
	%212 = fsub double -0.00000000e+00,  %211, !dbg !1713
	%213 = call double @llvm.fabs.f64 (double  %212), !dbg !1713
	store double  %213, double* %..inline.addr.550, align 8, !tbaa !1403, !dbg !1713
	store i32 0, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1714
	%214 = icmp sle i32  %200, 0, !dbg !1715
	br i1  %214, label %L..inline.10859, label %L.B0240, !dbg !1715
L.B0240:
	store i32  %200, i32* %.lcr014098.addr, align 4, !tbaa !1401, !dbg !1715
	br label %L..inline.10857
L..inline.10857:
	%215 = load i32, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1715
	%216 = sext i32  %215 to i64, !dbg !1715
	%217 = load double*, double** %..inline.addr.510, align 8, !tbaa !1398, !dbg !1715
	%218 = getelementptr double, double*  %217, i64  %216, !dbg !1715
	%219 = load double, double*  %218, align 8, !tbaa !1399, !dbg !1715
	%220 = fsub double -0.00000000e+00,  %219, !dbg !1715
	%221 = call double @llvm.fabs.f64 (double  %220), !dbg !1715
	store double  %221, double* %..inline.addr.570, align 8, !tbaa !1403, !dbg !1715
	%222 = fcmp une double  %221,  0.00000000000000000E+0, !dbg !1716
	br i1  %222, label %L..inline.10861, label %L.B0241, !dbg !1716
L.B0241:
	%223 = load double*, double** %..inline.addr.520, align 8, !tbaa !1398, !dbg !1717
	%224 = getelementptr double, double*  %223, i64  %216, !dbg !1717
	%225 = load double, double*  %224, align 8, !tbaa !1399, !dbg !1717
	call void @llvm.dbg.declare (metadata double* %result_rat.addr, metadata !1718, metadata !1366), !dbg !1617
	store double  %225, double* %result_rat.addr, align 8, !tbaa !1403, !dbg !1717
	call void @llvm.dbg.declare (metadata double* %dresult_rat.addr, metadata !1720, metadata !1366), !dbg !1617
	store double  0.00000000000000000E+0, double* %dresult_rat.addr, align 8, !tbaa !1403, !dbg !1719
	br label %L..inline.10916, !dbg !1721
L..inline.10861:
	%226 = load double, double* %..inline.addr.570, align 8, !tbaa !1403, !dbg !1722
	%227 = load double, double* %..inline.addr.550, align 8, !tbaa !1403, !dbg !1722
	%228 = fcmp uge double  %226,  %227, !dbg !1722
	br i1  %228, label %L..inline.10866, label %L.B0242, !dbg !1722
L.B0242:
	%229 = load i32, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1722
	store i32  %229, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !1722
	store double  %226, double* %..inline.addr.550, align 8, !tbaa !1403, !dbg !1723
	br label %L..inline.10866
L..inline.10866:
	%230 = load i32, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1724
	%231 = sext i32  %230 to i64, !dbg !1724
	%232 = load double*, double** %..inline.addr.520, align 8, !tbaa !1398, !dbg !1724
	%233 = getelementptr double, double*  %232, i64  %231, !dbg !1724
	%234 = load double, double*  %233, align 8, !tbaa !1399, !dbg !1724
	%235 = bitcast [20 x double]* %..inline.addr.600 to double*, !dbg !1724
	%236 = getelementptr double, double*  %235, i64  %231, !dbg !1724
	store double  %234, double*  %236, align 8, !tbaa !1399, !dbg !1724
	%237 = load double, double*  %233, align 8, !tbaa !1399, !dbg !1724
	%238 = fadd double  %237,  1.00000000000000008E-30, !dbg !1724
	%239 = bitcast [20 x double]* %..inline.addr.610 to double*, !dbg !1724
	%240 = getelementptr double, double*  %239, i64  %231, !dbg !1724
	store double  %238, double*  %240, align 8, !tbaa !1399, !dbg !1724

	%241 = add i32  %230, 1, !dbg !1726
	store i32  %241, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1726
	%242 = load i32, i32* %.lcr014098.addr, align 4, !tbaa !1401, !dbg !1726
	%243 = icmp slt i32  %241,  %242, !dbg !1726
	br i1  %243, label %L..inline.10857, label %L..inline.10859, !dbg !1726
L..inline.10859:
	%244 = load i32, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !1726

	%245 = sub i32  %244, 1, !dbg !1727
	store i32  %245, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !1727
	%246 = sext i32  %244 to i64, !dbg !1727
	%247 = load double*, double** %..inline.addr.520, align 8, !tbaa !1398, !dbg !1727
	%248 = getelementptr double, double*  %247, i64  %246, !dbg !1727
	%249 = load double, double*  %248, align 8, !tbaa !1399, !dbg !1727
	store double  %249, double* %result_rat.addr, align 8, !tbaa !1403, !dbg !1727
	store i32 1, i32* %..inline.addr.630, align 4, !tbaa !1401, !dbg !1727
	%250 = load i32, i32* %..inline.addr.530, align 4, !tbaa !1401, !dbg !1728
	%251 = icmp sge i32 1,  %250, !dbg !1728
	br i1  %251, label %L..inline.10871, label %L.B0243, !dbg !1728
L.B0243:
	%252 = sub i32  %250, 1, !dbg !1728
	store i32  %252, i32* %.x4098.addr, align 4, !tbaa !1401, !dbg !1728
	br label %L..inline.10870
L..inline.10870:
	store i32 0, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1728
	%253 = load i32, i32* %..inline.addr.530, align 4, !tbaa !1401, !dbg !1703
	%254 = load i32, i32* %..inline.addr.630, align 4, !tbaa !1401, !dbg !1703
	%255 = sub i32  %253,  %254, !dbg !1703
	%256 = icmp sle i32  %255, 0, !dbg !1703
	br i1  %256, label %L..inline.10873, label %L.B0244, !dbg !1703
L.B0244:
	store i32  %255, i32* %.lcr014099.addr, align 4, !tbaa !1401, !dbg !1703
	%257 = bitcast [20 x double]* %..inline.addr.600 to i8*, !dbg !1703
	%258 = getelementptr i8, i8*  %257, i64 8, !dbg !1703
	store i8*  %258, i8** %.G0037.addr, align 8, !tbaa !1398, !dbg !1703
	br label %L..inline.10872
L..inline.10872:
	%259 = load i8*, i8** %.G0037.addr, align 8, !tbaa !1398, !dbg !1703
	%260 = bitcast i8*  %259 to double*, !dbg !1703
	%261 = load double, double*  %260, align 8, !tbaa !1399, !dbg !1703
	%262 = load i32, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1703
	%263 = sext i32  %262 to i64, !dbg !1703
	%264 = bitcast [20 x double]* %..inline.addr.610 to double*, !dbg !1703
	%265 = getelementptr double, double*  %264, i64  %263, !dbg !1703
	%266 = load double, double*  %265, align 8, !tbaa !1399, !dbg !1703
	%267 = fsub double  %261,  %266, !dbg !1703
	store double  %267, double* %..inline.addr.650, align 8, !tbaa !1403, !dbg !1703
	%268 = load double*, double** %..inline.addr.510, align 8, !tbaa !1398, !dbg !1729
	%269 = getelementptr double, double*  %268, i64  %263, !dbg !1729
	%270 = load double, double*  %269, align 8, !tbaa !1399, !dbg !1729
	%271 = fmul double  %266,  %270, !dbg !1729
	%272 = load i32, i32* %..inline.addr.630, align 4, !tbaa !1401, !dbg !1729
	%273 = add i32  %262,  %272, !dbg !1729
	%274 = sext i32  %273 to i64, !dbg !1729
	%275 = getelementptr double, double*  %268, i64  %274, !dbg !1729
	%276 = load double, double*  %275, align 8, !tbaa !1399, !dbg !1729
	%277 = fdiv double  %271,  %276, !dbg !1729
	store double  %277, double* %..inline.addr.660, align 8, !tbaa !1403, !dbg !1729
	%278 = fsub double  %277,  %261, !dbg !1730
	store double  %278, double* %..inline.addr.670, align 8, !tbaa !1403, !dbg !1730
	%279 = fcmp une double  %278,  0.00000000000000000E+0, !dbg !1731
	br i1  %279, label %L..inline.10877, label %L.B0245, !dbg !1731
L.B0245:
	%280 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !1709
	%281 = icmp ne i8*  %280,  null, !dbg !1709
	br i1  %281, label %L..inline.10882, label %L.B0246, !dbg !1709
L.B0246:
	%282 = bitcast %struct._ZSo* @_ZSt4cerr to i8*, !dbg !1709
	%283 = bitcast %struct._ZSo* @_ZSt4cerr to i8**, !dbg !1709
	%284 = load i8*, i8**  %283, align 8, !tbaa !1398, !dbg !1709
	%285 = getelementptr i8, i8*  %284, i64 18446744073709551592, !dbg !1709
	%286 = bitcast i8*  %285 to i64*, !dbg !1709
	%287 = load i64, i64*  %286, align 8, !tbaa !1399, !dbg !1709
	%288 = getelementptr i8, i8*  %282, i64  %287, !dbg !1709
	%289 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.680 to i8**, !dbg !1709
	store i8*  %288, i8**  %289, align 8, !tbaa !1398, !dbg !1709
	%290 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.680, align 8, !tbaa !1398, !dbg !1709
	%291 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %290 to i8*, !dbg !1709
	%292 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.690 to i8**, !dbg !1709
	store i8*  %291, i8**  %292, align 8, !tbaa !1398, !dbg !1709
	%293 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.690, align 8, !tbaa !1398, !dbg !1709
	%294 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %293 to i8*, !dbg !1709
	%295 = getelementptr i8, i8*  %294, i64 32, !dbg !1709
	%296 = bitcast i8*  %295 to i32*, !dbg !1709
	%297 = load i32, i32*  %296, align 4, !tbaa !1399, !dbg !1709
	%298 = or i32  %297, 1, !dbg !1709
	call void  @_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate (%struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %290, i32  %298), !dbg !1709
	br label %L..inline.10903, !dbg !1709
L..inline.10882:
	%299 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !1709
	%300 = call i64  @strlen (i8*  %299) nounwind, !dbg !1709
	%301 = call %struct._ZSo*  @_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l (%struct._ZSo* @_ZSt4cerr, i8*  %299, i64  %300), !dbg !1709
	%302 = bitcast %struct._ZSo*  %301 to i8*, !dbg !1709
	%303 = bitcast %struct._ZSo** %.Q0010.addr to i8**, !dbg !1709
	store i8*  %302, i8**  %303, align 8, !tbaa !1398, !dbg !1709
	br label %L..inline.10903
L..inline.10903:
	call void  @exit (i32 111) nounwind noreturn, !dbg !1732
	br label %L..inline.10877
L..inline.10877:
	%304 = load double, double* %..inline.addr.650, align 8, !tbaa !1403, !dbg !1733
	%305 = load double, double* %..inline.addr.670, align 8, !tbaa !1403, !dbg !1733
	%306 = fdiv double  %304,  %305, !dbg !1733
	%307 = load i8*, i8** %.G0037.addr, align 8, !tbaa !1398, !dbg !1733
	%308 = bitcast i8*  %307 to double*, !dbg !1733
	%309 = load double, double*  %308, align 8, !tbaa !1399, !dbg !1733
	%310 = fmul double  %306,  %309, !dbg !1733
	%311 = load i32, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1733
	%312 = sext i32  %311 to i64, !dbg !1733
	%313 = bitcast [20 x double]* %..inline.addr.610 to double*, !dbg !1733
	%314 = getelementptr double, double*  %313, i64  %312, !dbg !1733
	store double  %310, double*  %314, align 8, !tbaa !1399, !dbg !1733
	%315 = load double, double* %..inline.addr.660, align 8, !tbaa !1403, !dbg !1734
	%316 = fmul double  %315,  %306, !dbg !1734
	%317 = getelementptr i8, i8*  %307, i64 18446744073709551608, !dbg !1734
	%318 = bitcast i8*  %317 to double*, !dbg !1734
	store double  %316, double*  %318, align 8, !tbaa !1399, !dbg !1734

	%319 = add i32  %311, 1, !dbg !1736
	store i32  %319, i32* %..inline.addr.560, align 4, !tbaa !1401, !dbg !1736
	%320 = getelementptr i8, i8*  %307, i64 8, !dbg !1703
	store i8*  %320, i8** %.G0037.addr, align 8, !tbaa !1398, !dbg !1703
	%321 = load i32, i32* %.lcr014099.addr, align 4, !tbaa !1401, !dbg !1736
	%322 = icmp slt i32  %319,  %321, !dbg !1736
	br i1  %322, label %L..inline.10872, label %L..inline.10873, !dbg !1736
L..inline.10873:
	%323 = load i32, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !1736
	%324 = add i32  %323, 1, !dbg !1736
	%325 = mul i32  %324, 2, !dbg !1736
	%326 = load i32, i32* %..inline.addr.530, align 4, !tbaa !1401, !dbg !1736
	%327 = load i32, i32* %..inline.addr.630, align 4, !tbaa !1401, !dbg !1736
	%328 = sub i32  %326,  %327, !dbg !1736
	%329 = icmp sge i32  %325,  %328, !dbg !1736
	br i1  %329, label %L..inline.10912, label %L.B0247, !dbg !1736
L.B0247:
	%330 = sext i32  %323 to i64, !dbg !1737
	%331 = bitcast [20 x double]* %..inline.addr.600 to i8*, !dbg !1737
	%332 = getelementptr i8, i8*  %331, i64 8, !dbg !1737
	%333 = bitcast i8*  %332 to double*, !dbg !1737
	%334 = getelementptr double, double*  %333, i64  %330, !dbg !1737
	%335 = load double, double*  %334, align 8, !tbaa !1399, !dbg !1737
	store double  %335, double* %..inline.addr.730, align 8, !tbaa !1403, !dbg !1737
	br label %L..inline.10914, !dbg !1737
L..inline.10912:
	%336 = load i32, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !1737

	%337 = sub i32  %336, 1, !dbg !1737
	store i32  %337, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !1737
	%338 = sext i32  %336 to i64, !dbg !1737
	%339 = bitcast [20 x double]* %..inline.addr.610 to double*, !dbg !1737
	%340 = getelementptr double, double*  %339, i64  %338, !dbg !1737
	%341 = load double, double*  %340, align 8, !tbaa !1399, !dbg !1737
	store double  %341, double* %..inline.addr.730, align 8, !tbaa !1403, !dbg !1737
	br label %L..inline.10914
L..inline.10914:
	%342 = load double, double* %..inline.addr.730, align 8, !tbaa !1403, !dbg !1737
	store double  %342, double* %dresult_rat.addr, align 8, !tbaa !1403, !dbg !1737
	%343 = load double, double* %result_rat.addr, align 8, !tbaa !1403, !dbg !1737
	%344 = fadd double  %342,  %343, !dbg !1737
	store double  %344, double* %result_rat.addr, align 8, !tbaa !1403, !dbg !1737
	%345 = load i32, i32* %..inline.addr.630, align 4, !tbaa !1401, !dbg !1737

	%346 = add i32  %345, 1, !dbg !1738
	store i32  %346, i32* %..inline.addr.630, align 4, !tbaa !1401, !dbg !1738
	%347 = load i32, i32* %.x4098.addr, align 4, !tbaa !1401, !dbg !1728
	%348 = sub i32  %347, 1, !dbg !1728
	store i32  %348, i32* %.x4098.addr, align 4, !tbaa !1401, !dbg !1728
	%349 = icmp sgt i32  %348, 0, !dbg !1738
	br i1  %349, label %L..inline.10870, label %L..inline.10871, !dbg !1738
L..inline.10871:
	br label %L..inline.10916
L..inline.10916:
	%350 = load double, double* %dresult_pol.addr, align 8, !tbaa !1403, !dbg !1739
	%351 = call double @llvm.fabs.f64 (double  %350), !dbg !1739
	store double  %351, double* %dresult_pol.addr, align 8, !tbaa !1403, !dbg !1739
	%352 = load double, double* %dresult_rat.addr, align 8, !tbaa !1403, !dbg !1740
	%353 = call double @llvm.fabs.f64 (double  %352), !dbg !1740
	store double  %353, double* %dresult_rat.addr, align 8, !tbaa !1403, !dbg !1740
	%354 = call double @llvm.maxnum.f64 (double  %351, double  %353), !dbg !1741
	%355 = fcmp uno double  %353,  %353, !dbg !1741
	%356 = select i1  %355, double  %351, double  %354, !dbg !1741
	store double  %356, double* %.I0003.addr, align 8, !tbaa !1403, !dbg !1741
	%357 = load double, double* %result_pol.addr, align 8, !tbaa !1403, !dbg !1742
	%358 = load double, double* %result_rat.addr, align 8, !tbaa !1403, !dbg !1742
	%359 = fsub double  %357,  %358, !dbg !1742
	%360 = call double @llvm.fabs.f64 (double  %359), !dbg !1742
	call void @llvm.dbg.declare (metadata double* %dresult_polrat.addr, metadata !1743, metadata !1366), !dbg !1617
	%361 = load double, double* %.I0003.addr, align 8, !tbaa !1403, !dbg !1744
	%362 = call double @llvm.maxnum.f64 (double  %361, double  %360), !dbg !1744
	%363 = fcmp uno double  %360,  %360, !dbg !1744
	%364 = select i1  %363, double  %361, double  %362, !dbg !1744
	call void @llvm.dbg.declare (metadata double* %dresult.addr, metadata !1745, metadata !1366), !dbg !1617
	%365 = load double, double* %absacc.addr, align 8, !tbaa !1403, !dbg !1746
	%366 = fcmp uge double  %364,  %365, !dbg !1746
	br i1  %366, label %L.B0047, label %L.B0248, !dbg !1746
L.B0248:
	%367 = load double, double* %result_pol.addr, align 8, !tbaa !1403, !dbg !1747
	%368 = load double, double* %result_rat.addr, align 8, !tbaa !1403, !dbg !1747
	%369 = fadd double  %367,  %368, !dbg !1747
	%370 = fmul double  %369,  5.00000000000000000E-1, !dbg !1747
	store double  %370, double* %result.addr, align 8, !tbaa !1403, !dbg !1747
	br label %L_T41511712_7633, !dbg !1748
L.B0047:
	br label %L.B0041
L.B0041:
	%371 = load i8*, i8** %.G0052.addr, align 8, !tbaa !1398, !dbg !1749
	%372 = getelementptr i8, i8*  %371, i64 18446744073709551608, !dbg !1749
	%373 = bitcast i8*  %372 to double*, !dbg !1749
	%374 = load double, double*  %373, align 8, !tbaa !1399, !dbg !1749
	%375 = bitcast i8*  %371 to double*, !dbg !1749
	store double  %374, double*  %375, align 8, !tbaa !1399, !dbg !1749
	%376 = load i8*, i8** %.G0050.addr, align 8, !tbaa !1398, !dbg !1750
	%377 = getelementptr i8, i8*  %376, i64 18446744073709551608, !dbg !1750
	%378 = bitcast i8*  %377 to double*, !dbg !1750
	%379 = load double, double*  %378, align 8, !tbaa !1399, !dbg !1750
	%380 = fmul double  %379,  2.50000000000000000E-1, !dbg !1750
	%381 = bitcast i8*  %376 to double*, !dbg !1750
	store double  %380, double*  %381, align 8, !tbaa !1399, !dbg !1750
	%382 = load i32, i32* %j.addr, align 4, !tbaa !1401, !dbg !1751

	%383 = add i32  %382, 1, !dbg !1751
	store i32  %383, i32* %j.addr, align 4, !tbaa !1401, !dbg !1751
	%384 = getelementptr i8, i8*  %371, i64 8, !dbg !1682
	store i8*  %384, i8** %.G0052.addr, align 8, !tbaa !1398, !dbg !1682
	%385 = getelementptr i8, i8*  %376, i64 8, !dbg !1682
	store i8*  %385, i8** %.G0050.addr, align 8, !tbaa !1398, !dbg !1682
	%386 = load i32, i32* %.x4097.addr, align 4, !tbaa !1401, !dbg !1682
	%387 = sub i32  %386, 1, !dbg !1682
	store i32  %387, i32* %.x4097.addr, align 4, !tbaa !1401, !dbg !1682
	%388 = icmp sgt i32  %387, 0, !dbg !1751
	br i1  %388, label %L.B0039, label %L.N0013, !dbg !1751
L.N0013:
	br label %L_T41511712_7633
L_T41511712_7633:
	%389 = load double, double* %result.addr, align 8, !tbaa !1403, !dbg !1752
	ret double  %389, !dbg !1752
L.N0011:
	br label %L.N0012
L.N0012:
	ret double 0.0
}

%struct.Tinty_f2D = type <{ %struct.T1DFunction, %struct.T2DFunction*, double, double, double}> 
%astruct.dt64 = type <{ i8*, i32}> 
%struct.T2DFunction = type <{ i32 (...)* (...)*}> 

define double @_Z7RombergRK11T2DFunctionddddd(%struct.T2DFunction* %func.arg, double %a.arg, double %b.arg, double %c.arg, double %d.arg, double %absacc.arg) #0 personality i8* bitcast (i32 (...)* @__gxx_personality_v0 to i8*) !dbg !1761 {
L.entry:
	%func.addr = alloca %struct.T2DFunction*, align 8
	%a.addr = alloca double, align 8
	%b.addr = alloca double, align 8
	%c.addr = alloca double, align 8
	%d.addr = alloca double, align 8
	%absacc.addr = alloca double, align 8
	%..inline.addr = alloca %struct.T2DFunction*, align 8
	%..inline.addr.1 = alloca double, align 8
	%..inline.addr.2 = alloca double, align 8
	%..inline.addr.3 = alloca double, align 8
	%_T40345360_7634.addr = alloca %struct.Tinty_f2D, align 8
	%..inline.addr.4 = alloca double, align 8
	%..inline.addr.5 = alloca double, align 8
	%..inline.addr.6 = alloca double, align 8
	%..inline.addr.7 = alloca i32, align 4
	%..inline.addr.8 = alloca [9 x double], align 8
	%..inline.addr.9 = alloca double, align 8
	%..inline.addr.230 = alloca i32, align 4
	%.x5121.addr = alloca i32, align 4
	%..inline.addr.250 = alloca [9 x double], align 8
	%.G0069.addr = alloca i8*, align 8
	%.G0067.addr = alloca i8*, align 8
	%..inline.addr.280 = alloca double, align 8
	%..inline.addr.290 = alloca double, align 8
	%..inline.addr.300 = alloca double*, align 8
	%.Q0011.addr = alloca double, align 8
	%.Q0012.addr = alloca double, align 8
	%..inline.addr.330 = alloca double, align 8
	%..inline.addr.340 = alloca double, align 8
	%..inline.addr.350 = alloca double, align 8
	%.x5122.addr = alloca i32, align 4
	%.Q0013.addr = alloca double, align 8
	%..inline.addr.380 = alloca double, align 8
	%..inline.addr.390 = alloca i32, align 4
	%..inline.addr.400 = alloca double*, align 8
	%..inline.addr.410 = alloca double*, align 8
	%..inline.addr.420 = alloca i32, align 4
	%..inline.addr.430 = alloca i32, align 4
	%..inline.addr.440 = alloca double, align 8
	%.ndi0014.addr = alloca i32, align 4
	%.lcr015122.addr = alloca i32, align 4
	%..inline.addr.470 = alloca double, align 8
	%..inline.addr.480 = alloca [20 x double], align 8
	%..inline.addr.490 = alloca [20 x double], align 8
	%..inline.addr.500 = alloca double, align 8
	%.TRP0006.addr = alloca i64, align 8
	%.ndk0028.addr = alloca i64, align 8
	%.lcr055122.addr = alloca i64, align 8
	%.TRP0007.addr = alloca i64, align 8
	%.ndk0030.addr = alloca i64, align 8
	%.lcr055123.addr = alloca i64, align 8
	%.r6.0563.addr = alloca i32, align 4
	%.G0054.addr = alloca i8*, align 8
	%..inline.addr.590 = alloca double, align 8
	%..inline.addr.600 = alloca double, align 8
	%..inline.addr.610 = alloca double, align 8
	%..inline.addr.620 = alloca double*, align 8
	%..inline.addr.630 = alloca double*, align 8
	%..inline.addr.640 = alloca i32, align 4
	%..inline.addr.650 = alloca i32, align 4
	%..inline.addr.660 = alloca double, align 8
	%..inline.addr.670 = alloca i32, align 4
	%..inline.addr.680 = alloca double, align 8
	%..inline.addr.690 = alloca double, align 8
	%..inline.addr.700 = alloca double, align 8
	%..inline.addr.710 = alloca [20 x double], align 8
	%..inline.addr.720 = alloca [20 x double], align 8
	%.D0028.addr = alloca i32, align 4
	%..inline.addr.740 = alloca i32, align 4
	%.lcr015123.addr = alloca i32, align 4
	%..inline.addr.760 = alloca double, align 8
	%..inline.addr.770 = alloca double, align 8
	%..inline.addr.780 = alloca double, align 8
	%..inline.addr.790 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.800 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.810 = alloca i32, align 4
	%..inline.addr.820 = alloca i64, align 8
	%.Q0014.addr = alloca %struct._ZSo*, align 8
	%..inline.addr.840 = alloca double, align 8
	%.D0029.addr = alloca i32, align 4
	%..inline.addr.860 = alloca double, align 8
	%..inline.addr.870 = alloca double, align 8
	%..inline.addr.880 = alloca double, align 8
	%_T40345656_7634.addr = alloca double, align 8
	%__caught_object_address.addr = alloca i8*, align 8
	%__catch_clause_number.addr = alloca i32, align 4

	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %func.addr, metadata !1853, metadata !1366), !dbg !1762
	store %struct.T2DFunction* %func.arg, %struct.T2DFunction** %func.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %func.addr, metadata !1854, metadata !1366), !dbg !1762
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1855, metadata !1366), !dbg !1762
	store double %a.arg, double* %a.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1856, metadata !1366), !dbg !1762
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1857, metadata !1366), !dbg !1762
	store double %b.arg, double* %b.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1858, metadata !1366), !dbg !1762
	call void @llvm.dbg.declare (metadata double* %c.addr, metadata !1859, metadata !1366), !dbg !1762
	store double %c.arg, double* %c.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %c.addr, metadata !1860, metadata !1366), !dbg !1762
	call void @llvm.dbg.declare (metadata double* %d.addr, metadata !1861, metadata !1366), !dbg !1762
	store double %d.arg, double* %d.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %d.addr, metadata !1862, metadata !1366), !dbg !1762
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !1863, metadata !1366), !dbg !1762
	store double %absacc.arg, double* %absacc.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !1864, metadata !1366), !dbg !1762
	%0 = load %struct.T2DFunction*, %struct.T2DFunction** %func.addr, align 8, !tbaa !1398, !dbg !1865
	%1 = bitcast %struct.T2DFunction*  %0 to i8*, !dbg !1865
	%2 = bitcast %struct.T2DFunction** %..inline.addr to i8**, !dbg !1865
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !1865
	%3 = load double, double* %c.addr, align 8, !tbaa !1403, !dbg !1865
	%4 = load double, double* %d.addr, align 8, !tbaa !1403, !dbg !1865
	%5 = load double, double* %absacc.addr, align 8, !tbaa !1403, !dbg !1865
	%6 = load double, double* %b.addr, align 8, !tbaa !1403, !dbg !1865
	%7 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !1865
	%8 = fsub double  %6,  %7, !dbg !1865
	%9 = fdiv double  %5,  %8, !dbg !1865
	%10 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !1866
	%11 = getelementptr i8, i8*  %10, i64 16, !dbg !1866
	%12 = bitcast %struct.Tinty_f2D* %_T40345360_7634.addr to i8**, !dbg !1866
	store i8*  %11, i8**  %12, align 8, !tbaa !1398, !dbg !1866
	%13 = bitcast [5 x i32 (...)* (...)*]* @_ZTV9Tinty_f2D to i8*, !dbg !1866
	%14 = getelementptr i8, i8*  %13, i64 16, !dbg !1866
	store i8*  %14, i8**  %12, align 8, !tbaa !1398, !dbg !1866
	%15 = load %struct.T2DFunction*, %struct.T2DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !1866
	%16 = bitcast %struct.T2DFunction*  %15 to i8*, !dbg !1866
	%17 = bitcast %struct.Tinty_f2D* %_T40345360_7634.addr to i8*, !dbg !1866
	%18 = getelementptr i8, i8*  %17, i64 8, !dbg !1866
	%19 = bitcast i8*  %18 to i8**, !dbg !1866
	store i8*  %16, i8**  %19, align 8, !tbaa !1398, !dbg !1866
	%20 = getelementptr i8, i8*  %17, i64 16, !dbg !1866
	%21 = bitcast i8*  %20 to double*, !dbg !1866
	store double  %3, double*  %21, align 8, !tbaa !1399, !dbg !1866
	%22 = getelementptr i8, i8*  %17, i64 24, !dbg !1866
	%23 = bitcast i8*  %22 to double*, !dbg !1866
	store double  %4, double*  %23, align 8, !tbaa !1399, !dbg !1866
	%24 = getelementptr i8, i8*  %17, i64 32, !dbg !1866
	%25 = bitcast i8*  %24 to double*, !dbg !1866
	store double  %9, double*  %25, align 8, !tbaa !1399, !dbg !1866
	store double  %7, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !1865
	store double  %6, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1865
	store double  %5, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1865
	store i32 0, i32* %..inline.addr.7, align 4, !tbaa !1401, !dbg !1865
	%26 = bitcast [9 x double]* %..inline.addr.8 to double*, !dbg !1874
	store double  1.00000000000000000E+0, double*  %26, align 8, !tbaa !1399, !dbg !1874
	store double  0.00000000000000000E+0, double* %..inline.addr.9, align 8, !tbaa !1403, !dbg !1875
	store i32 1, i32* %..inline.addr.230, align 4, !tbaa !1401, !dbg !1876
	store i32 8, i32* %.x5121.addr, align 4, !tbaa !1401, !dbg !1877
	%27 = bitcast [9 x double]* %..inline.addr.250 to i8*, !dbg !1877
	%28 = getelementptr i8, i8*  %27, i64 8, !dbg !1877
	store i8*  %28, i8** %.G0069.addr, align 8, !tbaa !1398, !dbg !1877
	%29 = bitcast [9 x double]* %..inline.addr.8 to i8*, !dbg !1877
	%30 = getelementptr i8, i8*  %29, i64 8, !dbg !1877
	store i8*  %30, i8** %.G0067.addr, align 8, !tbaa !1398, !dbg !1877
	br label %L..inline.11439
L..inline.11439:
	%31 = load double, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !1878
	store double  %31, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !1878
	%32 = load double, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1878
	store double  %32, double* %..inline.addr.290, align 8, !tbaa !1403, !dbg !1878
	%33 = bitcast [9 x double]* %..inline.addr.250 to i8*, !dbg !1878
	%34 = load i32, i32* %..inline.addr.230, align 4, !tbaa !1401, !dbg !1878
	%35 = sext i32  %34 to i64, !dbg !1878
	%36 = sub i64  %35, 1, !dbg !1878
	%37 = mul i64  %36, 8, !dbg !1878
	%38 = getelementptr i8, i8*  %33, i64  %37, !dbg !1878
	%39 = bitcast double** %..inline.addr.300 to i8**, !dbg !1878
	store i8*  %38, i8**  %39, align 8, !tbaa !1398, !dbg !1878
	%40 = icmp ne i32  %34, 1, !dbg !1878
	br i1  %40, label %L..inline.11453, label %L.B0272, !dbg !1878
L.B0272:
	%41 = bitcast %struct.Tinty_f2D* %_T40345360_7634.addr to %struct.T1DFunction*, !dbg !1878
	%42 = bitcast %struct.Tinty_f2D* %_T40345360_7634.addr to i8**, !dbg !1878
	%43 = load i8*, i8**  %42, align 8, !tbaa !1398, !dbg !1878
	%44 = bitcast i8*  %43 to double (%struct.T1DFunction*, double)**, !dbg !1878
	%45 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %44, align 8, !tbaa !1398, !dbg !1878
	%46 = invoke double  %45 (%struct.T1DFunction*  %41, double  %32)
		to label %L.B0273
		unwind label %L_T40346248_7634, !dbg !1878
L.B0273:
	store double  %46, double* %.Q0011.addr, align 8, !tbaa !1403, !dbg !1878
	%47 = load double, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !1878
	%48 = load i8*, i8**  %42, align 8, !tbaa !1398, !dbg !1878
	%49 = bitcast i8*  %48 to double (%struct.T1DFunction*, double)**, !dbg !1878
	%50 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %49, align 8, !tbaa !1398, !dbg !1878
	%51 = invoke double  %50 (%struct.T1DFunction*  %41, double  %47)
		to label %L.B0274
		unwind label %L_T40346248_7634, !dbg !1878
L.B0274:
	store double  %51, double* %.Q0012.addr, align 8, !tbaa !1403, !dbg !1878
	%52 = load double, double* %.Q0011.addr, align 8, !tbaa !1403, !dbg !1878
	%53 = fadd double  %52,  %51, !dbg !1878
	%54 = load double, double* %..inline.addr.290, align 8, !tbaa !1403, !dbg !1878
	%55 = load double, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !1878
	%56 = fsub double  %54,  %55, !dbg !1878
	%57 = fmul double  %56,  5.00000000000000000E-1, !dbg !1878
	%58 = fmul double  %53,  %57, !dbg !1878
	%59 = load double*, double** %..inline.addr.300, align 8, !tbaa !1398, !dbg !1878
	store double  %58, double*  %59, align 8, !tbaa !1399, !dbg !1878
	store i32 1, i32* %..inline.addr.7, align 4, !tbaa !1401, !dbg !1878
	br label %L..inline.11455, !dbg !1878
L..inline.11453:
	%60 = load double, double* %..inline.addr.5, align 8, !tbaa !1403, !dbg !1878
	%61 = load double, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !1878
	%62 = fsub double  %60,  %61, !dbg !1878
	%63 = load i32, i32* %..inline.addr.7, align 4, !tbaa !1401, !dbg !1878
	%64 = sitofp i32  %63 to double, !dbg !1878
	%65 = fdiv double  %62,  %64, !dbg !1878
	store double  %65, double* %..inline.addr.330, align 8, !tbaa !1403, !dbg !1878
	%66 = call double @llvm.fma.f64 (double  %65, double  5.00000000000000000E-1, double  %61), !dbg !1878
	store double  %66, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !1878
	store double  0.00000000000000000E+0, double* %..inline.addr.350, align 8, !tbaa !1403, !dbg !1878
	%67 = icmp sle i32  %63, 0, !dbg !1878
	br i1  %67, label %L..inline.11464, label %L.B0275, !dbg !1878
L.B0275:
	store i32  %63, i32* %.x5122.addr, align 4, !tbaa !1401, !dbg !1878
	br label %L..inline.11463
L..inline.11463:
	%68 = bitcast %struct.Tinty_f2D* %_T40345360_7634.addr to %struct.T1DFunction*, !dbg !1878
	%69 = load double, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !1878
	%70 = bitcast %struct.Tinty_f2D* %_T40345360_7634.addr to i8**, !dbg !1878
	%71 = load i8*, i8**  %70, align 8, !tbaa !1398, !dbg !1878
	%72 = bitcast i8*  %71 to double (%struct.T1DFunction*, double)**, !dbg !1878
	%73 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %72, align 8, !tbaa !1398, !dbg !1878
	%74 = invoke double  %73 (%struct.T1DFunction*  %68, double  %69)
		to label %L.B0276
		unwind label %L_T40346248_7634, !dbg !1878
L.B0276:
	store double  %74, double* %.Q0013.addr, align 8, !tbaa !1403, !dbg !1878
	%75 = load double, double* %..inline.addr.350, align 8, !tbaa !1403, !dbg !1878
	%76 = fadd double  %74,  %75, !dbg !1878
	store double  %76, double* %..inline.addr.350, align 8, !tbaa !1403, !dbg !1878
	%77 = load double, double* %..inline.addr.330, align 8, !tbaa !1403, !dbg !1878
	%78 = load double, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !1878
	%79 = fadd double  %77,  %78, !dbg !1878
	store double  %79, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !1878
	%80 = load i32, i32* %.x5122.addr, align 4, !tbaa !1401, !dbg !1878
	%81 = sub i32  %80, 1, !dbg !1878
	store i32  %81, i32* %.x5122.addr, align 4, !tbaa !1401, !dbg !1878
	%82 = icmp sgt i32  %81, 0, !dbg !1878
	br i1  %82, label %L..inline.11463, label %L..inline.11464, !dbg !1878
L..inline.11464:
	%83 = load double, double* %..inline.addr.290, align 8, !tbaa !1403, !dbg !1878
	%84 = load double, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !1878
	%85 = fsub double  %83,  %84, !dbg !1878
	%86 = load double, double* %..inline.addr.350, align 8, !tbaa !1403, !dbg !1878
	%87 = fmul double  %85,  %86, !dbg !1878
	%88 = load i32, i32* %..inline.addr.7, align 4, !tbaa !1401, !dbg !1878
	%89 = sitofp i32  %88 to double, !dbg !1878
	%90 = fdiv double  %87,  %89, !dbg !1878
	%91 = load double*, double** %..inline.addr.300, align 8, !tbaa !1398, !dbg !1878
	%92 = load double, double*  %91, align 8, !tbaa !1399, !dbg !1878
	%93 = fadd double  %90,  %92, !dbg !1878
	%94 = fmul double  %93,  5.00000000000000000E-1, !dbg !1878
	store double  %94, double*  %91, align 8, !tbaa !1399, !dbg !1878
	%95 = mul i32  %88, 2, !dbg !1878
	store i32  %95, i32* %..inline.addr.7, align 4, !tbaa !1401, !dbg !1878
	br label %L..inline.11455
L..inline.11455:
	%96 = load i8*, i8** %.G0069.addr, align 8, !tbaa !1398, !dbg !1878
	%97 = getelementptr i8, i8*  %96, i64 18446744073709551608, !dbg !1878
	%98 = bitcast i8*  %97 to double*, !dbg !1878
	%99 = load double, double*  %98, align 8, !tbaa !1399, !dbg !1878
	store double  %99, double* %..inline.addr.380, align 8, !tbaa !1403, !dbg !1878
	%100 = load i32, i32* %..inline.addr.230, align 4, !tbaa !1401, !dbg !1879
	%101 = icmp slt i32  %100, 2, !dbg !1879
	br i1  %101, label %L..inline.11466, label %L.B0277, !dbg !1879
L.B0277:
	%102 = icmp slt i32  %100, 5, !dbg !1880
	%103 = select i1  %102, i32  %100, i32 5, !dbg !1880
	store i32  %103, i32* %..inline.addr.390, align 4, !tbaa !1401, !dbg !1880
	%104 = bitcast [9 x double]* %..inline.addr.8 to i8*, !dbg !1881
	%105 = sub i32  %100,  %103, !dbg !1881
	%106 = sext i32  %105 to i64, !dbg !1881
	%107 = mul i64  %106, 8, !dbg !1881
	%108 = getelementptr i8, i8*  %104, i64  %107, !dbg !1881
	%109 = bitcast double** %..inline.addr.400 to i8**, !dbg !1881
	store i8*  %108, i8**  %109, align 8, !tbaa !1398, !dbg !1881
	%110 = bitcast [9 x double]* %..inline.addr.250 to i8*, !dbg !1881
	%111 = mul i64  %106, 8, !dbg !1881
	%112 = getelementptr i8, i8*  %110, i64  %111, !dbg !1881
	%113 = bitcast double** %..inline.addr.410 to i8**, !dbg !1881
	store i8*  %112, i8**  %113, align 8, !tbaa !1398, !dbg !1881
	store i32  %103, i32* %..inline.addr.420, align 4, !tbaa !1401, !dbg !1881
	store i32 0, i32* %..inline.addr.430, align 4, !tbaa !1401, !dbg !1881
	%114 = load double*, double** %..inline.addr.400, align 8, !tbaa !1398, !dbg !1881
	%115 = load double, double*  %114, align 8, !tbaa !1399, !dbg !1881
	%116 = fsub double -0.00000000e+00,  %115, !dbg !1881
	%117 = call double @llvm.fabs.f64 (double  %116), !dbg !1881
	store double  %117, double* %..inline.addr.440, align 8, !tbaa !1403, !dbg !1881
	store i32 0, i32* %.ndi0014.addr, align 4, !tbaa !1401, !dbg !1882
	%118 = icmp slt i32  %100, 5, !dbg !1882
	%119 = select i1  %118, i32  %100, i32 5, !dbg !1882
	%120 = icmp sle i32  %119, 0, !dbg !1882
	br i1  %120, label %L.B0267, label %L.B0278, !dbg !1882
L.B0278:
	%121 = icmp slt i32  %100, 5, !dbg !1881
	%122 = select i1  %121, i32  %100, i32 5, !dbg !1881
	store i32  %122, i32* %.lcr015122.addr, align 4, !tbaa !1401, !dbg !1881
	br label %L.B0266
L.B0266:
	%123 = load i32, i32* %.ndi0014.addr, align 4, !tbaa !1401, !dbg !1882
	%124 = sext i32  %123 to i64, !dbg !1882
	%125 = load double*, double** %..inline.addr.400, align 8, !tbaa !1398, !dbg !1882
	%126 = getelementptr double, double*  %125, i64  %124, !dbg !1882
	%127 = load double, double*  %126, align 8, !tbaa !1399, !dbg !1882
	%128 = fsub double -0.00000000e+00,  %127, !dbg !1882
	%129 = call double @llvm.fabs.f64 (double  %128), !dbg !1882
	store double  %129, double* %..inline.addr.470, align 8, !tbaa !1403, !dbg !1882
	%130 = load double, double* %..inline.addr.440, align 8, !tbaa !1403, !dbg !1882
	%131 = fcmp uge double  %129,  %130, !dbg !1882
	br i1  %131, label %L..inline.11488, label %L.B0279, !dbg !1882
L.B0279:
	store i32  %123, i32* %..inline.addr.430, align 4, !tbaa !1401, !dbg !1882
	store double  %129, double* %..inline.addr.440, align 8, !tbaa !1403, !dbg !1882
	br label %L..inline.11488
L..inline.11488:
	%132 = load i32, i32* %.ndi0014.addr, align 4, !tbaa !1401, !dbg !1882
	%133 = sext i32  %132 to i64, !dbg !1882
	%134 = load double*, double** %..inline.addr.410, align 8, !tbaa !1398, !dbg !1882
	%135 = getelementptr double, double*  %134, i64  %133, !dbg !1882
	%136 = load double, double*  %135, align 8, !tbaa !1399, !dbg !1882
	%137 = bitcast [20 x double]* %..inline.addr.480 to double*, !dbg !1882
	%138 = getelementptr double, double*  %137, i64  %133, !dbg !1882
	store double  %136, double*  %138, align 8, !tbaa !1399, !dbg !1882
	%139 = bitcast [20 x double]* %..inline.addr.490 to double*, !dbg !1882
	%140 = getelementptr double, double*  %139, i64  %133, !dbg !1882
	store double  %136, double*  %140, align 8, !tbaa !1399, !dbg !1882
	%141 = add i32  %132, 1, !dbg !1882
	store i32  %141, i32* %.ndi0014.addr, align 4, !tbaa !1401, !dbg !1882
	%142 = load i32, i32* %.lcr015122.addr, align 4, !tbaa !1401, !dbg !1882
	%143 = icmp slt i32  %141,  %142, !dbg !1882
	br i1  %143, label %L.B0266, label %L.B0267, !dbg !1882
L.B0267:
	%144 = load i32, i32* %..inline.addr.430, align 4, !tbaa !1401, !dbg !1881
	%145 = sext i32  %144 to i64, !dbg !1881
	%146 = load double*, double** %..inline.addr.410, align 8, !tbaa !1398, !dbg !1881
	%147 = getelementptr double, double*  %146, i64  %145, !dbg !1881
	%148 = load double, double*  %147, align 8, !tbaa !1399, !dbg !1881
	store double  %148, double* %..inline.addr.500, align 8, !tbaa !1403, !dbg !1881

	%149 = sub i32  %144, 1, !dbg !1881
	store i32  %149, i32* %..inline.addr.430, align 4, !tbaa !1401, !dbg !1881
	%150 = load i32, i32* %..inline.addr.420, align 4, !tbaa !1401, !dbg !1865
	%151 = sext i32  %150 to i64, !dbg !1865
	%152 = sub i64  %151, 1, !dbg !1865
	store i64  %152, i64* %.TRP0006.addr, align 8, !tbaa !1449, !dbg !1865
	store i64 0, i64* %.ndk0028.addr, align 8, !tbaa !1451, !dbg !1882
	%153 = icmp sle i64  %152, 0, !dbg !1882
	br i1  %153, label %L.B0269, label %L.B0280, !dbg !1882
L.B0280:
	store i64  %152, i64* %.lcr055122.addr, align 8, !tbaa !1451, !dbg !1881
	br label %L.B0268
L.B0268:
	%154 = load i32, i32* %..inline.addr.420, align 4, !tbaa !1401, !dbg !1882
	%155 = load i64, i64* %.ndk0028.addr, align 8, !tbaa !1451, !dbg !1882
	%156 = trunc i64  %155 to i32, !dbg !1882
	%157 = sub i32  %154,  %156, !dbg !1882
	%158 = sext i32  %157 to i64, !dbg !1882
	%159 = sub i64  %158, 1, !dbg !1882
	store i64  %159, i64* %.TRP0007.addr, align 8, !tbaa !1449, !dbg !1882
	store i64 0, i64* %.ndk0030.addr, align 8, !tbaa !1451, !dbg !1882
	%160 = icmp sle i64  %159, 0, !dbg !1882
	br i1  %160, label %L.B0271, label %L.B0281, !dbg !1882
L.B0281:
	store i64  %159, i64* %.lcr055123.addr, align 8, !tbaa !1451, !dbg !1881
	store i32  %156, i32* %.r6.0563.addr, align 4, !tbaa !1401, !dbg !1881
	%161 = bitcast [20 x double]* %..inline.addr.490 to i8*, !dbg !1881
	store i8*  %161, i8** %.G0054.addr, align 8, !tbaa !1398, !dbg !1881
	br label %L.B0270
L.B0270:
	%162 = load i64, i64* %.ndk0030.addr, align 8, !tbaa !1451, !dbg !1882
	%163 = trunc i64  %162 to i32, !dbg !1882
	%164 = sext i32  %163 to i64, !dbg !1882
	%165 = load double*, double** %..inline.addr.400, align 8, !tbaa !1398, !dbg !1882
	%166 = getelementptr double, double*  %165, i64  %164, !dbg !1882
	%167 = load double, double*  %166, align 8, !tbaa !1399, !dbg !1882
	%168 = load i8*, i8** %.G0054.addr, align 8, !tbaa !1398, !dbg !1882
	%169 = getelementptr i8, i8*  %168, i64 8, !dbg !1882
	%170 = bitcast i8*  %169 to double*, !dbg !1882
	%171 = load double, double*  %170, align 8, !tbaa !1399, !dbg !1882
	%172 = bitcast [20 x double]* %..inline.addr.480 to double*, !dbg !1882
	%173 = getelementptr double, double*  %172, i64  %164, !dbg !1882
	%174 = load double, double*  %173, align 8, !tbaa !1399, !dbg !1882
	%175 = fsub double  %171,  %174, !dbg !1882
	%176 = load i32, i32* %.r6.0563.addr, align 4, !tbaa !1401, !dbg !1882
	%177 = add i32  %176,  %163, !dbg !1882
	%178 = sext i32  %177 to i64, !dbg !1882
	%179 = bitcast double*  %165 to i8*, !dbg !1882
	%180 = getelementptr i8, i8*  %179, i64 8, !dbg !1882
	%181 = bitcast i8*  %180 to double*, !dbg !1882
	%182 = getelementptr double, double*  %181, i64  %178, !dbg !1882
	%183 = load double, double*  %182, align 8, !tbaa !1399, !dbg !1882
	%184 = fsub double  %167,  %183, !dbg !1882
	%185 = fdiv double  %175,  %184, !dbg !1882
	%186 = fmul double  %185,  %183, !dbg !1882
	store double  %186, double*  %173, align 8, !tbaa !1399, !dbg !1882
	%187 = fmul double  %167,  %185, !dbg !1882
	%188 = bitcast i8*  %168 to double*, !dbg !1882
	store double  %187, double*  %188, align 8, !tbaa !1399, !dbg !1882
	%189 = add i64  %162, 1, !dbg !1882
	store i64  %189, i64* %.ndk0030.addr, align 8, !tbaa !1451, !dbg !1882
	store i8*  %169, i8** %.G0054.addr, align 8, !tbaa !1398, !dbg !1881
	%190 = load i64, i64* %.lcr055123.addr, align 8, !tbaa !1451, !dbg !1882
	%191 = icmp slt i64  %189,  %190, !dbg !1882
	br i1  %191, label %L.B0270, label %L.B0271, !dbg !1882
L.B0271:
	%192 = load i32, i32* %..inline.addr.430, align 4, !tbaa !1401, !dbg !1882
	%193 = add i32  %192, 1, !dbg !1882
	%194 = mul i32  %193, 2, !dbg !1882
	%195 = load i32, i32* %..inline.addr.420, align 4, !tbaa !1401, !dbg !1882
	%196 = load i64, i64* %.ndk0028.addr, align 8, !tbaa !1451, !dbg !1882
	%197 = trunc i64  %196 to i32, !dbg !1882
	%198 = sub i32  %195,  %197, !dbg !1882
	%199 = sub i32  %198, 1, !dbg !1882
	%200 = icmp sge i32  %194,  %199, !dbg !1882
	br i1  %200, label %L..inline.11507, label %L.B0282, !dbg !1882
L.B0282:
	%201 = sext i32  %192 to i64, !dbg !1882
	%202 = bitcast [20 x double]* %..inline.addr.490 to i8*, !dbg !1882
	%203 = getelementptr i8, i8*  %202, i64 8, !dbg !1882
	%204 = bitcast i8*  %203 to double*, !dbg !1882
	%205 = getelementptr double, double*  %204, i64  %201, !dbg !1882
	%206 = load double, double*  %205, align 8, !tbaa !1399, !dbg !1882
	store double  %206, double* %..inline.addr.9, align 8, !tbaa !1403, !dbg !1882
	br label %L..inline.11509, !dbg !1882
L..inline.11507:
	%207 = load i32, i32* %..inline.addr.430, align 4, !tbaa !1401, !dbg !1882
	%208 = sext i32  %207 to i64, !dbg !1882
	%209 = bitcast [20 x double]* %..inline.addr.480 to double*, !dbg !1882
	%210 = getelementptr double, double*  %209, i64  %208, !dbg !1882
	%211 = load double, double*  %210, align 8, !tbaa !1399, !dbg !1882
	store double  %211, double* %..inline.addr.9, align 8, !tbaa !1403, !dbg !1882

	%212 = sub i32  %207, 1, !dbg !1882
	store i32  %212, i32* %..inline.addr.430, align 4, !tbaa !1401, !dbg !1882
	br label %L..inline.11509
L..inline.11509:
	%213 = load double, double* %..inline.addr.9, align 8, !tbaa !1403, !dbg !1882
	%214 = load double, double* %..inline.addr.500, align 8, !tbaa !1403, !dbg !1882
	%215 = fadd double  %213,  %214, !dbg !1882
	store double  %215, double* %..inline.addr.500, align 8, !tbaa !1403, !dbg !1882
	%216 = load i64, i64* %.ndk0028.addr, align 8, !tbaa !1451, !dbg !1882
	%217 = add i64  %216, 1, !dbg !1882
	store i64  %217, i64* %.ndk0028.addr, align 8, !tbaa !1451, !dbg !1882
	%218 = load i64, i64* %.lcr055122.addr, align 8, !tbaa !1451, !dbg !1882
	%219 = icmp slt i64  %217,  %218, !dbg !1882
	br i1  %219, label %L.B0268, label %L.B0269, !dbg !1882
L.B0269:
	%220 = bitcast [9 x double]* %..inline.addr.8 to i8*, !dbg !1882
	%221 = load i32, i32* %..inline.addr.230, align 4, !tbaa !1401, !dbg !1882
	%222 = load i32, i32* %..inline.addr.390, align 4, !tbaa !1401, !dbg !1882
	%223 = sub i32  %221,  %222, !dbg !1882
	%224 = sext i32  %223 to i64, !dbg !1882
	%225 = mul i64  %224, 8, !dbg !1882
	%226 = getelementptr i8, i8*  %220, i64  %225, !dbg !1882
	%227 = bitcast double** %..inline.addr.620 to i8**, !dbg !1882
	store i8*  %226, i8**  %227, align 8, !tbaa !1398, !dbg !1882
	%228 = bitcast [9 x double]* %..inline.addr.250 to i8*, !dbg !1882
	%229 = mul i64  %224, 8, !dbg !1882
	%230 = getelementptr i8, i8*  %228, i64  %229, !dbg !1882
	%231 = bitcast double** %..inline.addr.630 to i8**, !dbg !1882
	store i8*  %230, i8**  %231, align 8, !tbaa !1398, !dbg !1882
	store i32  %222, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !1882
	store i32 0, i32* %..inline.addr.650, align 4, !tbaa !1401, !dbg !1882
	%232 = load double*, double** %..inline.addr.620, align 8, !tbaa !1398, !dbg !1882
	%233 = load double, double*  %232, align 8, !tbaa !1399, !dbg !1882
	%234 = fsub double -0.00000000e+00,  %233, !dbg !1882
	%235 = call double @llvm.fabs.f64 (double  %234), !dbg !1882
	store double  %235, double* %..inline.addr.660, align 8, !tbaa !1403, !dbg !1882
	store i32 0, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	%236 = icmp sle i32  %222, 0, !dbg !1882
	br i1  %236, label %L..inline.11525, label %L.B0283, !dbg !1882
L.B0283:
	store i32  %222, i32* %.lcr015122.addr, align 4, !tbaa !1401, !dbg !1882
	br label %L..inline.11524
L..inline.11524:
	%237 = load i32, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	%238 = sext i32  %237 to i64, !dbg !1882
	%239 = load double*, double** %..inline.addr.620, align 8, !tbaa !1398, !dbg !1882
	%240 = getelementptr double, double*  %239, i64  %238, !dbg !1882
	%241 = load double, double*  %240, align 8, !tbaa !1399, !dbg !1882
	%242 = fsub double -0.00000000e+00,  %241, !dbg !1882
	%243 = call double @llvm.fabs.f64 (double  %242), !dbg !1882
	store double  %243, double* %..inline.addr.680, align 8, !tbaa !1403, !dbg !1882
	%244 = fcmp une double  %243,  0.00000000000000000E+0, !dbg !1882
	br i1  %244, label %L..inline.11527, label %L.B0284, !dbg !1882
L.B0284:
	%245 = load double*, double** %..inline.addr.630, align 8, !tbaa !1398, !dbg !1882
	%246 = getelementptr double, double*  %245, i64  %238, !dbg !1882
	%247 = load double, double*  %246, align 8, !tbaa !1399, !dbg !1882
	store double  %247, double* %..inline.addr.690, align 8, !tbaa !1403, !dbg !1882
	store double  0.00000000000000000E+0, double* %..inline.addr.700, align 8, !tbaa !1403, !dbg !1882
	br label %L..inline.11530, !dbg !1882
L..inline.11527:
	%248 = load double, double* %..inline.addr.680, align 8, !tbaa !1403, !dbg !1882
	%249 = load double, double* %..inline.addr.660, align 8, !tbaa !1403, !dbg !1882
	%250 = fcmp uge double  %248,  %249, !dbg !1882
	br i1  %250, label %L..inline.11532, label %L.B0285, !dbg !1882
L.B0285:
	%251 = load i32, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	store i32  %251, i32* %..inline.addr.650, align 4, !tbaa !1401, !dbg !1882
	store double  %248, double* %..inline.addr.660, align 8, !tbaa !1403, !dbg !1882
	br label %L..inline.11532
L..inline.11532:
	%252 = load i32, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	%253 = sext i32  %252 to i64, !dbg !1882
	%254 = load double*, double** %..inline.addr.630, align 8, !tbaa !1398, !dbg !1882
	%255 = getelementptr double, double*  %254, i64  %253, !dbg !1882
	%256 = load double, double*  %255, align 8, !tbaa !1399, !dbg !1882
	%257 = bitcast [20 x double]* %..inline.addr.710 to double*, !dbg !1882
	%258 = getelementptr double, double*  %257, i64  %253, !dbg !1882
	store double  %256, double*  %258, align 8, !tbaa !1399, !dbg !1882
	%259 = load double, double*  %255, align 8, !tbaa !1399, !dbg !1882
	%260 = fadd double  %259,  1.00000000000000008E-30, !dbg !1882
	%261 = bitcast [20 x double]* %..inline.addr.720 to double*, !dbg !1882
	%262 = getelementptr double, double*  %261, i64  %253, !dbg !1882
	store double  %260, double*  %262, align 8, !tbaa !1399, !dbg !1882

	%263 = add i32  %252, 1, !dbg !1882
	store i32  %263, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	%264 = load i32, i32* %.lcr015122.addr, align 4, !tbaa !1401, !dbg !1882
	%265 = icmp slt i32  %263,  %264, !dbg !1882
	br i1  %265, label %L..inline.11524, label %L..inline.11525, !dbg !1882
L..inline.11525:
	%266 = load i32, i32* %..inline.addr.650, align 4, !tbaa !1401, !dbg !1882

	%267 = sub i32  %266, 1, !dbg !1882
	store i32  %267, i32* %..inline.addr.650, align 4, !tbaa !1401, !dbg !1882
	%268 = sext i32  %266 to i64, !dbg !1882
	%269 = load double*, double** %..inline.addr.630, align 8, !tbaa !1398, !dbg !1882
	%270 = getelementptr double, double*  %269, i64  %268, !dbg !1882
	%271 = load double, double*  %270, align 8, !tbaa !1399, !dbg !1882
	store double  %271, double* %..inline.addr.690, align 8, !tbaa !1403, !dbg !1882
	store i32 1, i32* %..inline.addr.740, align 4, !tbaa !1401, !dbg !1882
	%272 = load i32, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !1882
	%273 = icmp sge i32 1,  %272, !dbg !1882
	br i1  %273, label %L..inline.11537, label %L.B0286, !dbg !1882
L.B0286:
	%274 = sub i32  %272, 1, !dbg !1882
	store i32  %274, i32* %.x5122.addr, align 4, !tbaa !1401, !dbg !1882
	br label %L..inline.11536
L..inline.11536:
	store i32 0, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	%275 = load i32, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !1882
	%276 = load i32, i32* %..inline.addr.740, align 4, !tbaa !1401, !dbg !1882
	%277 = sub i32  %275,  %276, !dbg !1882
	%278 = icmp sle i32  %277, 0, !dbg !1882
	br i1  %278, label %L..inline.11539, label %L.B0287, !dbg !1882
L.B0287:
	store i32  %277, i32* %.lcr015123.addr, align 4, !tbaa !1401, !dbg !1882
	%279 = bitcast [20 x double]* %..inline.addr.710 to i8*, !dbg !1882
	%280 = getelementptr i8, i8*  %279, i64 8, !dbg !1882
	store i8*  %280, i8** %.G0054.addr, align 8, !tbaa !1398, !dbg !1882
	br label %L..inline.11538
L..inline.11538:
	%281 = load i8*, i8** %.G0054.addr, align 8, !tbaa !1398, !dbg !1882
	%282 = bitcast i8*  %281 to double*, !dbg !1882
	%283 = load double, double*  %282, align 8, !tbaa !1399, !dbg !1882
	%284 = load i32, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	%285 = sext i32  %284 to i64, !dbg !1882
	%286 = bitcast [20 x double]* %..inline.addr.720 to double*, !dbg !1882
	%287 = getelementptr double, double*  %286, i64  %285, !dbg !1882
	%288 = load double, double*  %287, align 8, !tbaa !1399, !dbg !1882
	%289 = fsub double  %283,  %288, !dbg !1882
	store double  %289, double* %..inline.addr.760, align 8, !tbaa !1403, !dbg !1882
	%290 = load double*, double** %..inline.addr.620, align 8, !tbaa !1398, !dbg !1882
	%291 = getelementptr double, double*  %290, i64  %285, !dbg !1882
	%292 = load double, double*  %291, align 8, !tbaa !1399, !dbg !1882
	%293 = fmul double  %288,  %292, !dbg !1882
	%294 = load i32, i32* %..inline.addr.740, align 4, !tbaa !1401, !dbg !1882
	%295 = add i32  %284,  %294, !dbg !1882
	%296 = sext i32  %295 to i64, !dbg !1882
	%297 = getelementptr double, double*  %290, i64  %296, !dbg !1882
	%298 = load double, double*  %297, align 8, !tbaa !1399, !dbg !1882
	%299 = fdiv double  %293,  %298, !dbg !1882
	store double  %299, double* %..inline.addr.770, align 8, !tbaa !1403, !dbg !1882
	%300 = fsub double  %299,  %283, !dbg !1882
	store double  %300, double* %..inline.addr.780, align 8, !tbaa !1403, !dbg !1882
	%301 = fcmp une double  %300,  0.00000000000000000E+0, !dbg !1882
	br i1  %301, label %L..inline.11543, label %L.B0288, !dbg !1882
L.B0288:
	%302 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !1882
	%303 = icmp ne i8*  %302,  null, !dbg !1882
	br i1  %303, label %L..inline.11548, label %L.B0289, !dbg !1882
L.B0289:
	%304 = bitcast %struct._ZSo* @_ZSt4cerr to i8*, !dbg !1882
	%305 = bitcast %struct._ZSo* @_ZSt4cerr to i8**, !dbg !1882
	%306 = load i8*, i8**  %305, align 8, !tbaa !1398, !dbg !1882
	%307 = getelementptr i8, i8*  %306, i64 18446744073709551592, !dbg !1882
	%308 = bitcast i8*  %307 to i64*, !dbg !1882
	%309 = load i64, i64*  %308, align 8, !tbaa !1399, !dbg !1882
	%310 = getelementptr i8, i8*  %304, i64  %309, !dbg !1882
	%311 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.790 to i8**, !dbg !1882
	store i8*  %310, i8**  %311, align 8, !tbaa !1398, !dbg !1882
	%312 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.790, align 8, !tbaa !1398, !dbg !1882
	%313 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %312 to i8*, !dbg !1882
	%314 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.800 to i8**, !dbg !1882
	store i8*  %313, i8**  %314, align 8, !tbaa !1398, !dbg !1882
	%315 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.800, align 8, !tbaa !1398, !dbg !1882
	%316 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %315 to i8*, !dbg !1882
	%317 = getelementptr i8, i8*  %316, i64 32, !dbg !1882
	%318 = bitcast i8*  %317 to i32*, !dbg !1882
	%319 = load i32, i32*  %318, align 4, !tbaa !1399, !dbg !1882
	%320 = or i32  %319, 1, !dbg !1882
	invoke void  @_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate (%struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %312, i32  %320)
		to label %L.B0290
		unwind label %L_T40346248_7634, !dbg !1882
L.B0290:
	br label %L..inline.11569, !dbg !1882
L..inline.11548:
	%321 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !1882
	%322 = call i64  @strlen (i8*  %321) nounwind, !dbg !1882
	%323 = invoke %struct._ZSo*  @_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l (%struct._ZSo* @_ZSt4cerr, i8*  %321, i64  %322)
		to label %L.B0291
		unwind label %L_T40346248_7634, !dbg !1882
L.B0291:
	%324 = bitcast %struct._ZSo*  %323 to i8*, !dbg !1882
	%325 = bitcast %struct._ZSo** %.Q0014.addr to i8**, !dbg !1882
	store i8*  %324, i8**  %325, align 8, !tbaa !1398, !dbg !1882
	br label %L..inline.11569
L..inline.11569:
	call void  @exit (i32 111) nounwind noreturn, !dbg !1882
	br label %L..inline.11543
L..inline.11543:
	%326 = load double, double* %..inline.addr.760, align 8, !tbaa !1403, !dbg !1882
	%327 = load double, double* %..inline.addr.780, align 8, !tbaa !1403, !dbg !1882
	%328 = fdiv double  %326,  %327, !dbg !1882
	%329 = load i8*, i8** %.G0054.addr, align 8, !tbaa !1398, !dbg !1882
	%330 = bitcast i8*  %329 to double*, !dbg !1882
	%331 = load double, double*  %330, align 8, !tbaa !1399, !dbg !1882
	%332 = fmul double  %328,  %331, !dbg !1882
	%333 = load i32, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	%334 = sext i32  %333 to i64, !dbg !1882
	%335 = bitcast [20 x double]* %..inline.addr.720 to double*, !dbg !1882
	%336 = getelementptr double, double*  %335, i64  %334, !dbg !1882
	store double  %332, double*  %336, align 8, !tbaa !1399, !dbg !1882
	%337 = load double, double* %..inline.addr.770, align 8, !tbaa !1403, !dbg !1882
	%338 = fmul double  %337,  %328, !dbg !1882
	%339 = getelementptr i8, i8*  %329, i64 18446744073709551608, !dbg !1882
	%340 = bitcast i8*  %339 to double*, !dbg !1882
	store double  %338, double*  %340, align 8, !tbaa !1399, !dbg !1882

	%341 = add i32  %333, 1, !dbg !1882
	store i32  %341, i32* %..inline.addr.670, align 4, !tbaa !1401, !dbg !1882
	%342 = getelementptr i8, i8*  %329, i64 8, !dbg !1882
	store i8*  %342, i8** %.G0054.addr, align 8, !tbaa !1398, !dbg !1882
	%343 = load i32, i32* %.lcr015123.addr, align 4, !tbaa !1401, !dbg !1882
	%344 = icmp slt i32  %341,  %343, !dbg !1882
	br i1  %344, label %L..inline.11538, label %L..inline.11539, !dbg !1882
L..inline.11539:
	%345 = load i32, i32* %..inline.addr.650, align 4, !tbaa !1401, !dbg !1882
	%346 = add i32  %345, 1, !dbg !1882
	%347 = mul i32  %346, 2, !dbg !1882
	%348 = load i32, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !1882
	%349 = load i32, i32* %..inline.addr.740, align 4, !tbaa !1401, !dbg !1882
	%350 = sub i32  %348,  %349, !dbg !1882
	%351 = icmp sge i32  %347,  %350, !dbg !1882
	br i1  %351, label %L..inline.11578, label %L.B0292, !dbg !1882
L.B0292:
	%352 = sext i32  %345 to i64, !dbg !1882
	%353 = bitcast [20 x double]* %..inline.addr.710 to i8*, !dbg !1882
	%354 = getelementptr i8, i8*  %353, i64 8, !dbg !1882
	%355 = bitcast i8*  %354 to double*, !dbg !1882
	%356 = getelementptr double, double*  %355, i64  %352, !dbg !1882
	%357 = load double, double*  %356, align 8, !tbaa !1399, !dbg !1882
	store double  %357, double* %..inline.addr.840, align 8, !tbaa !1403, !dbg !1882
	br label %L..inline.11580, !dbg !1882
L..inline.11578:
	%358 = load i32, i32* %..inline.addr.650, align 4, !tbaa !1401, !dbg !1882

	%359 = sub i32  %358, 1, !dbg !1882
	store i32  %359, i32* %..inline.addr.650, align 4, !tbaa !1401, !dbg !1882
	%360 = sext i32  %358 to i64, !dbg !1882
	%361 = bitcast [20 x double]* %..inline.addr.720 to double*, !dbg !1882
	%362 = getelementptr double, double*  %361, i64  %360, !dbg !1882
	%363 = load double, double*  %362, align 8, !tbaa !1399, !dbg !1882
	store double  %363, double* %..inline.addr.840, align 8, !tbaa !1403, !dbg !1882
	br label %L..inline.11580
L..inline.11580:
	%364 = load double, double* %..inline.addr.840, align 8, !tbaa !1403, !dbg !1882
	store double  %364, double* %..inline.addr.700, align 8, !tbaa !1403, !dbg !1882
	%365 = load double, double* %..inline.addr.690, align 8, !tbaa !1403, !dbg !1882
	%366 = fadd double  %364,  %365, !dbg !1882
	store double  %366, double* %..inline.addr.690, align 8, !tbaa !1403, !dbg !1882
	%367 = load i32, i32* %..inline.addr.740, align 4, !tbaa !1401, !dbg !1882

	%368 = add i32  %367, 1, !dbg !1882
	store i32  %368, i32* %..inline.addr.740, align 4, !tbaa !1401, !dbg !1882
	%369 = load i32, i32* %.x5122.addr, align 4, !tbaa !1401, !dbg !1882
	%370 = sub i32  %369, 1, !dbg !1882
	store i32  %370, i32* %.x5122.addr, align 4, !tbaa !1401, !dbg !1882
	%371 = icmp sgt i32  %370, 0, !dbg !1882
	br i1  %371, label %L..inline.11536, label %L..inline.11537, !dbg !1882
L..inline.11537:
	br label %L..inline.11530
L..inline.11530:
	%372 = load double, double* %..inline.addr.9, align 8, !tbaa !1403, !dbg !1882
	%373 = call double @llvm.fabs.f64 (double  %372), !dbg !1882
	store double  %373, double* %..inline.addr.9, align 8, !tbaa !1403, !dbg !1882
	%374 = load double, double* %..inline.addr.700, align 8, !tbaa !1403, !dbg !1883
	%375 = call double @llvm.fabs.f64 (double  %374), !dbg !1883
	store double  %375, double* %..inline.addr.700, align 8, !tbaa !1403, !dbg !1883
	%376 = call double @llvm.maxnum.f64 (double  %373, double  %375), !dbg !1884
	%377 = fcmp uno double  %375,  %375, !dbg !1884
	%378 = select i1  %377, double  %373, double  %376, !dbg !1884
	%379 = load double, double* %..inline.addr.500, align 8, !tbaa !1403, !dbg !1885
	%380 = load double, double* %..inline.addr.690, align 8, !tbaa !1403, !dbg !1885
	%381 = fsub double  %379,  %380, !dbg !1885
	%382 = call double @llvm.fabs.f64 (double  %381), !dbg !1885
	%383 = call double @llvm.maxnum.f64 (double  %378, double  %382), !dbg !1886
	%384 = fcmp uno double  %382,  %382, !dbg !1886
	%385 = select i1  %384, double  %378, double  %383, !dbg !1886
	%386 = load double, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !1887
	%387 = fcmp uge double  %385,  %386, !dbg !1887
	br i1  %387, label %L..inline.11588, label %L.B0293, !dbg !1887
L.B0293:
	%388 = fadd double  %379,  %380, !dbg !1887
	%389 = fmul double  %388,  5.00000000000000000E-1, !dbg !1887
	store double  %389, double* %..inline.addr.380, align 8, !tbaa !1403, !dbg !1887
	br label %L..inline.11589, !dbg !1888
L..inline.11588:
	br label %L..inline.11466
L..inline.11466:
	%390 = load i8*, i8** %.G0069.addr, align 8, !tbaa !1398, !dbg !1889
	%391 = getelementptr i8, i8*  %390, i64 18446744073709551608, !dbg !1889
	%392 = bitcast i8*  %391 to double*, !dbg !1889
	%393 = load double, double*  %392, align 8, !tbaa !1399, !dbg !1889
	%394 = bitcast i8*  %390 to double*, !dbg !1889
	store double  %393, double*  %394, align 8, !tbaa !1399, !dbg !1889
	%395 = load i8*, i8** %.G0067.addr, align 8, !tbaa !1398, !dbg !1889
	%396 = getelementptr i8, i8*  %395, i64 18446744073709551608, !dbg !1889
	%397 = bitcast i8*  %396 to double*, !dbg !1889
	%398 = load double, double*  %397, align 8, !tbaa !1399, !dbg !1889
	%399 = fmul double  %398,  2.50000000000000000E-1, !dbg !1889
	%400 = bitcast i8*  %395 to double*, !dbg !1889
	store double  %399, double*  %400, align 8, !tbaa !1399, !dbg !1889
	%401 = load i32, i32* %..inline.addr.230, align 4, !tbaa !1401, !dbg !1890

	%402 = add i32  %401, 1, !dbg !1891
	store i32  %402, i32* %..inline.addr.230, align 4, !tbaa !1401, !dbg !1891
	%403 = getelementptr i8, i8*  %390, i64 8, !dbg !1877
	store i8*  %403, i8** %.G0069.addr, align 8, !tbaa !1398, !dbg !1877
	%404 = getelementptr i8, i8*  %395, i64 8, !dbg !1877
	store i8*  %404, i8** %.G0067.addr, align 8, !tbaa !1398, !dbg !1877
	%405 = load i32, i32* %.x5121.addr, align 4, !tbaa !1401, !dbg !1877
	%406 = sub i32  %405, 1, !dbg !1877
	store i32  %406, i32* %.x5121.addr, align 4, !tbaa !1401, !dbg !1877
	%407 = icmp sgt i32  %406, 0, !dbg !1891
	br i1  %407, label %L..inline.11439, label %L.N0018, !dbg !1891
L.N0018:
	br label %L..inline.11589
L..inline.11589:
	%408 = load double, double* %..inline.addr.380, align 8, !tbaa !1403, !dbg !1892
	%409 = bitcast [5 x i32 (...)* (...)*]* @_ZTV9Tinty_f2D to i8*, !dbg !1865
	%410 = getelementptr i8, i8*  %409, i64 16, !dbg !1865
	%411 = bitcast %struct.Tinty_f2D* %_T40345360_7634.addr to i8**, !dbg !1865
	store i8*  %410, i8**  %411, align 8, !tbaa !1398, !dbg !1865
	%412 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !1865
	%413 = getelementptr i8, i8*  %412, i64 16, !dbg !1865
	store i8*  %413, i8**  %411, align 8, !tbaa !1398, !dbg !1865
	ret double  %408, !dbg !1865
L.B0294:
	br label %L.R0005, !dbg !1865
L_T40346248_7634:
	%414 = landingpad %astruct.dt64
	cleanup
	%415 = extractvalue %astruct.dt64  %414, 0, !dbg !1865
	store i8*  %415, i8** %__caught_object_address.addr, align 1, !tbaa !1398, !dbg !1865
	%416 = extractvalue %astruct.dt64  %414, 1, !dbg !1865
	store i32  %416, i32* %__catch_clause_number.addr, align 1, !tbaa !1399, !dbg !1865
	%417 = bitcast [5 x i32 (...)* (...)*]* @_ZTV9Tinty_f2D to i8*, !dbg !1865
	%418 = getelementptr i8, i8*  %417, i64 16, !dbg !1865
	%419 = bitcast %struct.Tinty_f2D* %_T40345360_7634.addr to i8**, !dbg !1865
	store i8*  %418, i8**  %419, align 8, !tbaa !1398, !dbg !1865
	%420 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !1865
	%421 = getelementptr i8, i8*  %420, i64 16, !dbg !1865
	store i8*  %421, i8**  %419, align 8, !tbaa !1398, !dbg !1865
	%422 = load i32, i32* %__catch_clause_number.addr, align 4, !tbaa !1399, !dbg !1865
	%423 = load i8*, i8** %__caught_object_address.addr, align 8, !tbaa !1398, !dbg !1865
	%424 = insertvalue %astruct.dt64 undef, i8*  %423, 0, !dbg !1865
	%425 = insertvalue %astruct.dt64  %424, i32  %422, 1, !dbg !1865
	resume %astruct.dt64  %425 , !dbg !1865
L.R0005:
	ret double 0.0
}

%struct.T3DFunction = type <{ i32 (...)* (...)*}> 
%struct.Tintxy_f3D = type <{ %struct.T1DFunction, %struct.T3DFunction*, double, double, double, double, double}> 

define double @_Z7RombergRK11T3DFunctionddddddd(%struct.T3DFunction* %func.arg, double %a.arg, double %b.arg, double %c.arg, double %d.arg, double %e.arg, double %f.arg, double %absacc.arg) #0 personality i8* bitcast (i32 (...)* @__gxx_personality_v0 to i8*) !dbg !1901 {
L.entry:
	%func.addr = alloca %struct.T3DFunction*, align 8
	%a.addr = alloca double, align 8
	%b.addr = alloca double, align 8
	%c.addr = alloca double, align 8
	%d.addr = alloca double, align 8
	%e.addr = alloca double, align 8
	%f.addr = alloca double, align 8
	%absacc.addr = alloca double, align 8
	%..inline.addr = alloca %struct.T3DFunction*, align 8
	%..inline.addr.1 = alloca double, align 8
	%..inline.addr.2 = alloca double, align 8
	%..inline.addr.3 = alloca double, align 8
	%..inline.addr.4 = alloca double, align 8
	%..inline.addr.5 = alloca double, align 8
	%_T40345952_7635.addr = alloca %struct.Tintxy_f3D, align 8
	%..inline.addr.6 = alloca double, align 8
	%..inline.addr.7 = alloca double, align 8
	%..inline.addr.8 = alloca double, align 8
	%..inline.addr.9 = alloca i32, align 4
	%..inline.addr.270 = alloca [9 x double], align 8
	%..inline.addr.280 = alloca double, align 8
	%..inline.addr.290 = alloca i32, align 4
	%.x6145.addr = alloca i32, align 4
	%..inline.addr.310 = alloca [9 x double], align 8
	%.G0086.addr = alloca i8*, align 8
	%.G0084.addr = alloca i8*, align 8
	%..inline.addr.340 = alloca double, align 8
	%..inline.addr.350 = alloca double, align 8
	%..inline.addr.360 = alloca double*, align 8
	%.Q0015.addr = alloca double, align 8
	%.Q0016.addr = alloca double, align 8
	%..inline.addr.390 = alloca double, align 8
	%..inline.addr.400 = alloca double, align 8
	%..inline.addr.410 = alloca double, align 8
	%.x6146.addr = alloca i32, align 4
	%.Q0017.addr = alloca double, align 8
	%..inline.addr.440 = alloca double, align 8
	%..inline.addr.450 = alloca i32, align 4
	%..inline.addr.460 = alloca double*, align 8
	%..inline.addr.470 = alloca double*, align 8
	%..inline.addr.480 = alloca i32, align 4
	%..inline.addr.490 = alloca i32, align 4
	%..inline.addr.500 = alloca double, align 8
	%.ndi0018.addr = alloca i32, align 4
	%.lcr016146.addr = alloca i32, align 4
	%..inline.addr.530 = alloca double, align 8
	%..inline.addr.540 = alloca [20 x double], align 8
	%..inline.addr.550 = alloca [20 x double], align 8
	%..inline.addr.560 = alloca double, align 8
	%.TRP0008.addr = alloca i64, align 8
	%.ndk0036.addr = alloca i64, align 8
	%.lcr056146.addr = alloca i64, align 8
	%.TRP0009.addr = alloca i64, align 8
	%.ndk0038.addr = alloca i64, align 8
	%.lcr056147.addr = alloca i64, align 8
	%.r7.0577.addr = alloca i32, align 4
	%.G0071.addr = alloca i8*, align 8
	%..inline.addr.650 = alloca double, align 8
	%..inline.addr.660 = alloca double, align 8
	%..inline.addr.670 = alloca double, align 8
	%..inline.addr.680 = alloca double*, align 8
	%..inline.addr.690 = alloca double*, align 8
	%..inline.addr.700 = alloca i32, align 4
	%..inline.addr.710 = alloca i32, align 4
	%..inline.addr.720 = alloca double, align 8
	%..inline.addr.730 = alloca i32, align 4
	%..inline.addr.740 = alloca double, align 8
	%..inline.addr.750 = alloca double, align 8
	%..inline.addr.760 = alloca double, align 8
	%..inline.addr.770 = alloca [20 x double], align 8
	%..inline.addr.780 = alloca [20 x double], align 8
	%.D0037.addr = alloca i32, align 4
	%..inline.addr.800 = alloca i32, align 4
	%.lcr016147.addr = alloca i32, align 4
	%..inline.addr.820 = alloca double, align 8
	%..inline.addr.830 = alloca double, align 8
	%..inline.addr.840 = alloca double, align 8
	%..inline.addr.850 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.860 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.870 = alloca i32, align 4
	%..inline.addr.880 = alloca i64, align 8
	%.Q0018.addr = alloca %struct._ZSo*, align 8
	%..inline.addr.900 = alloca double, align 8
	%.D0038.addr = alloca i32, align 4
	%..inline.addr.920 = alloca double, align 8
	%..inline.addr.930 = alloca double, align 8
	%..inline.addr.940 = alloca double, align 8
	%_T40346248_7635.addr = alloca double, align 8
	%__caught_object_address.addr = alloca i8*, align 8
	%__catch_clause_number.addr = alloca i32, align 4

	call void @llvm.dbg.declare (metadata %struct.T3DFunction** %func.addr, metadata !1993, metadata !1366), !dbg !1902
	store %struct.T3DFunction* %func.arg, %struct.T3DFunction** %func.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T3DFunction** %func.addr, metadata !1994, metadata !1366), !dbg !1902
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1995, metadata !1366), !dbg !1902
	store double %a.arg, double* %a.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %a.addr, metadata !1996, metadata !1366), !dbg !1902
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1997, metadata !1366), !dbg !1902
	store double %b.arg, double* %b.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %b.addr, metadata !1998, metadata !1366), !dbg !1902
	call void @llvm.dbg.declare (metadata double* %c.addr, metadata !1999, metadata !1366), !dbg !1902
	store double %c.arg, double* %c.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %c.addr, metadata !2000, metadata !1366), !dbg !1902
	call void @llvm.dbg.declare (metadata double* %d.addr, metadata !2001, metadata !1366), !dbg !1902
	store double %d.arg, double* %d.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %d.addr, metadata !2002, metadata !1366), !dbg !1902
	call void @llvm.dbg.declare (metadata double* %e.addr, metadata !2003, metadata !1366), !dbg !1902
	store double %e.arg, double* %e.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %e.addr, metadata !2004, metadata !1366), !dbg !1902
	call void @llvm.dbg.declare (metadata double* %f.addr, metadata !2005, metadata !1366), !dbg !1902
	store double %f.arg, double* %f.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %f.addr, metadata !2006, metadata !1366), !dbg !1902
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !2007, metadata !1366), !dbg !1902
	store double %absacc.arg, double* %absacc.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !2008, metadata !1366), !dbg !1902
	%0 = load %struct.T3DFunction*, %struct.T3DFunction** %func.addr, align 8, !tbaa !1398, !dbg !2009
	%1 = bitcast %struct.T3DFunction*  %0 to i8*, !dbg !2009
	%2 = bitcast %struct.T3DFunction** %..inline.addr to i8**, !dbg !2009
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2009
	%3 = load double, double* %a.addr, align 8, !tbaa !1403, !dbg !2009
	%4 = load double, double* %b.addr, align 8, !tbaa !1403, !dbg !2009
	%5 = load double, double* %c.addr, align 8, !tbaa !1403, !dbg !2009
	%6 = load double, double* %d.addr, align 8, !tbaa !1403, !dbg !2009
	%7 = load double, double* %absacc.addr, align 8, !tbaa !1403, !dbg !2009
	%8 = load double, double* %f.addr, align 8, !tbaa !1403, !dbg !2009
	%9 = load double, double* %e.addr, align 8, !tbaa !1403, !dbg !2009
	%10 = fsub double  %8,  %9, !dbg !2009
	%11 = fdiv double  %7,  %10, !dbg !2009
	%12 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2010
	%13 = getelementptr i8, i8*  %12, i64 16, !dbg !2010
	%14 = bitcast %struct.Tintxy_f3D* %_T40345952_7635.addr to i8**, !dbg !2010
	store i8*  %13, i8**  %14, align 8, !tbaa !1398, !dbg !2010
	%15 = bitcast [5 x i32 (...)* (...)*]* @_ZTV10Tintxy_f3D to i8*, !dbg !2010
	%16 = getelementptr i8, i8*  %15, i64 16, !dbg !2010
	store i8*  %16, i8**  %14, align 8, !tbaa !1398, !dbg !2010
	%17 = load %struct.T3DFunction*, %struct.T3DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !2010
	%18 = bitcast %struct.T3DFunction*  %17 to i8*, !dbg !2010
	%19 = bitcast %struct.Tintxy_f3D* %_T40345952_7635.addr to i8*, !dbg !2010
	%20 = getelementptr i8, i8*  %19, i64 8, !dbg !2010
	%21 = bitcast i8*  %20 to i8**, !dbg !2010
	store i8*  %18, i8**  %21, align 8, !tbaa !1398, !dbg !2010
	%22 = getelementptr i8, i8*  %19, i64 16, !dbg !2010
	%23 = bitcast i8*  %22 to double*, !dbg !2010
	store double  %3, double*  %23, align 8, !tbaa !1399, !dbg !2010
	%24 = getelementptr i8, i8*  %19, i64 24, !dbg !2010
	%25 = bitcast i8*  %24 to double*, !dbg !2010
	store double  %4, double*  %25, align 8, !tbaa !1399, !dbg !2010
	%26 = getelementptr i8, i8*  %19, i64 32, !dbg !2010
	%27 = bitcast i8*  %26 to double*, !dbg !2010
	store double  %5, double*  %27, align 8, !tbaa !1399, !dbg !2010
	%28 = getelementptr i8, i8*  %19, i64 40, !dbg !2010
	%29 = bitcast i8*  %28 to double*, !dbg !2010
	store double  %6, double*  %29, align 8, !tbaa !1399, !dbg !2010
	%30 = getelementptr i8, i8*  %19, i64 48, !dbg !2010
	%31 = bitcast i8*  %30 to double*, !dbg !2010
	store double  %11, double*  %31, align 8, !tbaa !1399, !dbg !2010
	store double  %9, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !2009
	store double  %8, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2009
	store double  %7, double* %..inline.addr.8, align 8, !tbaa !1403, !dbg !2009
	store i32 0, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !2009
	%32 = bitcast [9 x double]* %..inline.addr.270 to double*, !dbg !2013
	store double  1.00000000000000000E+0, double*  %32, align 8, !tbaa !1399, !dbg !2013
	store double  0.00000000000000000E+0, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2014
	store i32 1, i32* %..inline.addr.290, align 4, !tbaa !1401, !dbg !2015
	store i32 8, i32* %.x6145.addr, align 4, !tbaa !1401, !dbg !2016
	%33 = bitcast [9 x double]* %..inline.addr.310 to i8*, !dbg !2016
	%34 = getelementptr i8, i8*  %33, i64 8, !dbg !2016
	store i8*  %34, i8** %.G0086.addr, align 8, !tbaa !1398, !dbg !2016
	%35 = bitcast [9 x double]* %..inline.addr.270 to i8*, !dbg !2016
	%36 = getelementptr i8, i8*  %35, i64 8, !dbg !2016
	store i8*  %36, i8** %.G0084.addr, align 8, !tbaa !1398, !dbg !2016
	br label %L..inline.12178
L..inline.12178:
	%37 = load double, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !2017
	store double  %37, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !2017
	%38 = load double, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2017
	store double  %38, double* %..inline.addr.350, align 8, !tbaa !1403, !dbg !2017
	%39 = bitcast [9 x double]* %..inline.addr.310 to i8*, !dbg !2017
	%40 = load i32, i32* %..inline.addr.290, align 4, !tbaa !1401, !dbg !2017
	%41 = sext i32  %40 to i64, !dbg !2017
	%42 = sub i64  %41, 1, !dbg !2017
	%43 = mul i64  %42, 8, !dbg !2017
	%44 = getelementptr i8, i8*  %39, i64  %43, !dbg !2017
	%45 = bitcast double** %..inline.addr.360 to i8**, !dbg !2017
	store i8*  %44, i8**  %45, align 8, !tbaa !1398, !dbg !2017
	%46 = icmp ne i32  %40, 1, !dbg !2017
	br i1  %46, label %L..inline.12192, label %L.B0318, !dbg !2017
L.B0318:
	%47 = bitcast %struct.Tintxy_f3D* %_T40345952_7635.addr to %struct.T1DFunction*, !dbg !2017
	%48 = bitcast %struct.Tintxy_f3D* %_T40345952_7635.addr to i8**, !dbg !2017
	%49 = load i8*, i8**  %48, align 8, !tbaa !1398, !dbg !2017
	%50 = bitcast i8*  %49 to double (%struct.T1DFunction*, double)**, !dbg !2017
	%51 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %50, align 8, !tbaa !1398, !dbg !2017
	%52 = invoke double  %51 (%struct.T1DFunction*  %47, double  %38)
		to label %L.B0319
		unwind label %L_T41505256_7635, !dbg !2017
L.B0319:
	store double  %52, double* %.Q0015.addr, align 8, !tbaa !1403, !dbg !2017
	%53 = load double, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !2017
	%54 = load i8*, i8**  %48, align 8, !tbaa !1398, !dbg !2017
	%55 = bitcast i8*  %54 to double (%struct.T1DFunction*, double)**, !dbg !2017
	%56 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %55, align 8, !tbaa !1398, !dbg !2017
	%57 = invoke double  %56 (%struct.T1DFunction*  %47, double  %53)
		to label %L.B0320
		unwind label %L_T41505256_7635, !dbg !2017
L.B0320:
	store double  %57, double* %.Q0016.addr, align 8, !tbaa !1403, !dbg !2017
	%58 = load double, double* %.Q0015.addr, align 8, !tbaa !1403, !dbg !2017
	%59 = fadd double  %58,  %57, !dbg !2017
	%60 = load double, double* %..inline.addr.350, align 8, !tbaa !1403, !dbg !2017
	%61 = load double, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !2017
	%62 = fsub double  %60,  %61, !dbg !2017
	%63 = fmul double  %62,  5.00000000000000000E-1, !dbg !2017
	%64 = fmul double  %59,  %63, !dbg !2017
	%65 = load double*, double** %..inline.addr.360, align 8, !tbaa !1398, !dbg !2017
	store double  %64, double*  %65, align 8, !tbaa !1399, !dbg !2017
	store i32 1, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !2017
	br label %L..inline.12194, !dbg !2017
L..inline.12192:
	%66 = load double, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2017
	%67 = load double, double* %..inline.addr.6, align 8, !tbaa !1403, !dbg !2017
	%68 = fsub double  %66,  %67, !dbg !2017
	%69 = load i32, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !2017
	%70 = sitofp i32  %69 to double, !dbg !2017
	%71 = fdiv double  %68,  %70, !dbg !2017
	store double  %71, double* %..inline.addr.390, align 8, !tbaa !1403, !dbg !2017
	%72 = call double @llvm.fma.f64 (double  %71, double  5.00000000000000000E-1, double  %67), !dbg !2017
	store double  %72, double* %..inline.addr.400, align 8, !tbaa !1403, !dbg !2017
	store double  0.00000000000000000E+0, double* %..inline.addr.410, align 8, !tbaa !1403, !dbg !2017
	%73 = icmp sle i32  %69, 0, !dbg !2017
	br i1  %73, label %L..inline.12203, label %L.B0321, !dbg !2017
L.B0321:
	store i32  %69, i32* %.x6146.addr, align 4, !tbaa !1401, !dbg !2017
	br label %L..inline.12202
L..inline.12202:
	%74 = bitcast %struct.Tintxy_f3D* %_T40345952_7635.addr to %struct.T1DFunction*, !dbg !2017
	%75 = load double, double* %..inline.addr.400, align 8, !tbaa !1403, !dbg !2017
	%76 = bitcast %struct.Tintxy_f3D* %_T40345952_7635.addr to i8**, !dbg !2017
	%77 = load i8*, i8**  %76, align 8, !tbaa !1398, !dbg !2017
	%78 = bitcast i8*  %77 to double (%struct.T1DFunction*, double)**, !dbg !2017
	%79 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %78, align 8, !tbaa !1398, !dbg !2017
	%80 = invoke double  %79 (%struct.T1DFunction*  %74, double  %75)
		to label %L.B0322
		unwind label %L_T41505256_7635, !dbg !2017
L.B0322:
	store double  %80, double* %.Q0017.addr, align 8, !tbaa !1403, !dbg !2017
	%81 = load double, double* %..inline.addr.410, align 8, !tbaa !1403, !dbg !2017
	%82 = fadd double  %80,  %81, !dbg !2017
	store double  %82, double* %..inline.addr.410, align 8, !tbaa !1403, !dbg !2017
	%83 = load double, double* %..inline.addr.390, align 8, !tbaa !1403, !dbg !2017
	%84 = load double, double* %..inline.addr.400, align 8, !tbaa !1403, !dbg !2017
	%85 = fadd double  %83,  %84, !dbg !2017
	store double  %85, double* %..inline.addr.400, align 8, !tbaa !1403, !dbg !2017
	%86 = load i32, i32* %.x6146.addr, align 4, !tbaa !1401, !dbg !2017
	%87 = sub i32  %86, 1, !dbg !2017
	store i32  %87, i32* %.x6146.addr, align 4, !tbaa !1401, !dbg !2017
	%88 = icmp sgt i32  %87, 0, !dbg !2017
	br i1  %88, label %L..inline.12202, label %L..inline.12203, !dbg !2017
L..inline.12203:
	%89 = load double, double* %..inline.addr.350, align 8, !tbaa !1403, !dbg !2017
	%90 = load double, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !2017
	%91 = fsub double  %89,  %90, !dbg !2017
	%92 = load double, double* %..inline.addr.410, align 8, !tbaa !1403, !dbg !2017
	%93 = fmul double  %91,  %92, !dbg !2017
	%94 = load i32, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !2017
	%95 = sitofp i32  %94 to double, !dbg !2017
	%96 = fdiv double  %93,  %95, !dbg !2017
	%97 = load double*, double** %..inline.addr.360, align 8, !tbaa !1398, !dbg !2017
	%98 = load double, double*  %97, align 8, !tbaa !1399, !dbg !2017
	%99 = fadd double  %96,  %98, !dbg !2017
	%100 = fmul double  %99,  5.00000000000000000E-1, !dbg !2017
	store double  %100, double*  %97, align 8, !tbaa !1399, !dbg !2017
	%101 = mul i32  %94, 2, !dbg !2017
	store i32  %101, i32* %..inline.addr.9, align 4, !tbaa !1401, !dbg !2017
	br label %L..inline.12194
L..inline.12194:
	%102 = load i8*, i8** %.G0086.addr, align 8, !tbaa !1398, !dbg !2017
	%103 = getelementptr i8, i8*  %102, i64 18446744073709551608, !dbg !2017
	%104 = bitcast i8*  %103 to double*, !dbg !2017
	%105 = load double, double*  %104, align 8, !tbaa !1399, !dbg !2017
	store double  %105, double* %..inline.addr.440, align 8, !tbaa !1403, !dbg !2017
	%106 = load i32, i32* %..inline.addr.290, align 4, !tbaa !1401, !dbg !2018
	%107 = icmp slt i32  %106, 2, !dbg !2018
	br i1  %107, label %L..inline.12205, label %L.B0323, !dbg !2018
L.B0323:
	%108 = icmp slt i32  %106, 5, !dbg !2019
	%109 = select i1  %108, i32  %106, i32 5, !dbg !2019
	store i32  %109, i32* %..inline.addr.450, align 4, !tbaa !1401, !dbg !2019
	%110 = bitcast [9 x double]* %..inline.addr.270 to i8*, !dbg !2020
	%111 = sub i32  %106,  %109, !dbg !2020
	%112 = sext i32  %111 to i64, !dbg !2020
	%113 = mul i64  %112, 8, !dbg !2020
	%114 = getelementptr i8, i8*  %110, i64  %113, !dbg !2020
	%115 = bitcast double** %..inline.addr.460 to i8**, !dbg !2020
	store i8*  %114, i8**  %115, align 8, !tbaa !1398, !dbg !2020
	%116 = bitcast [9 x double]* %..inline.addr.310 to i8*, !dbg !2020
	%117 = mul i64  %112, 8, !dbg !2020
	%118 = getelementptr i8, i8*  %116, i64  %117, !dbg !2020
	%119 = bitcast double** %..inline.addr.470 to i8**, !dbg !2020
	store i8*  %118, i8**  %119, align 8, !tbaa !1398, !dbg !2020
	store i32  %109, i32* %..inline.addr.480, align 4, !tbaa !1401, !dbg !2020
	store i32 0, i32* %..inline.addr.490, align 4, !tbaa !1401, !dbg !2020
	%120 = load double*, double** %..inline.addr.460, align 8, !tbaa !1398, !dbg !2020
	%121 = load double, double*  %120, align 8, !tbaa !1399, !dbg !2020
	%122 = fsub double -0.00000000e+00,  %121, !dbg !2020
	%123 = call double @llvm.fabs.f64 (double  %122), !dbg !2020
	store double  %123, double* %..inline.addr.500, align 8, !tbaa !1403, !dbg !2020
	store i32 0, i32* %.ndi0018.addr, align 4, !tbaa !1401, !dbg !2021
	%124 = icmp slt i32  %106, 5, !dbg !2021
	%125 = select i1  %124, i32  %106, i32 5, !dbg !2021
	%126 = icmp sle i32  %125, 0, !dbg !2021
	br i1  %126, label %L.B0313, label %L.B0324, !dbg !2021
L.B0324:
	%127 = icmp slt i32  %106, 5, !dbg !2020
	%128 = select i1  %127, i32  %106, i32 5, !dbg !2020
	store i32  %128, i32* %.lcr016146.addr, align 4, !tbaa !1401, !dbg !2020
	br label %L.B0312
L.B0312:
	%129 = load i32, i32* %.ndi0018.addr, align 4, !tbaa !1401, !dbg !2021
	%130 = sext i32  %129 to i64, !dbg !2021
	%131 = load double*, double** %..inline.addr.460, align 8, !tbaa !1398, !dbg !2021
	%132 = getelementptr double, double*  %131, i64  %130, !dbg !2021
	%133 = load double, double*  %132, align 8, !tbaa !1399, !dbg !2021
	%134 = fsub double -0.00000000e+00,  %133, !dbg !2021
	%135 = call double @llvm.fabs.f64 (double  %134), !dbg !2021
	store double  %135, double* %..inline.addr.530, align 8, !tbaa !1403, !dbg !2021
	%136 = load double, double* %..inline.addr.500, align 8, !tbaa !1403, !dbg !2021
	%137 = fcmp uge double  %135,  %136, !dbg !2021
	br i1  %137, label %L..inline.12227, label %L.B0325, !dbg !2021
L.B0325:
	store i32  %129, i32* %..inline.addr.490, align 4, !tbaa !1401, !dbg !2021
	store double  %135, double* %..inline.addr.500, align 8, !tbaa !1403, !dbg !2021
	br label %L..inline.12227
L..inline.12227:
	%138 = load i32, i32* %.ndi0018.addr, align 4, !tbaa !1401, !dbg !2021
	%139 = sext i32  %138 to i64, !dbg !2021
	%140 = load double*, double** %..inline.addr.470, align 8, !tbaa !1398, !dbg !2021
	%141 = getelementptr double, double*  %140, i64  %139, !dbg !2021
	%142 = load double, double*  %141, align 8, !tbaa !1399, !dbg !2021
	%143 = bitcast [20 x double]* %..inline.addr.540 to double*, !dbg !2021
	%144 = getelementptr double, double*  %143, i64  %139, !dbg !2021
	store double  %142, double*  %144, align 8, !tbaa !1399, !dbg !2021
	%145 = bitcast [20 x double]* %..inline.addr.550 to double*, !dbg !2021
	%146 = getelementptr double, double*  %145, i64  %139, !dbg !2021
	store double  %142, double*  %146, align 8, !tbaa !1399, !dbg !2021
	%147 = add i32  %138, 1, !dbg !2021
	store i32  %147, i32* %.ndi0018.addr, align 4, !tbaa !1401, !dbg !2021
	%148 = load i32, i32* %.lcr016146.addr, align 4, !tbaa !1401, !dbg !2021
	%149 = icmp slt i32  %147,  %148, !dbg !2021
	br i1  %149, label %L.B0312, label %L.B0313, !dbg !2021
L.B0313:
	%150 = load i32, i32* %..inline.addr.490, align 4, !tbaa !1401, !dbg !2020
	%151 = sext i32  %150 to i64, !dbg !2020
	%152 = load double*, double** %..inline.addr.470, align 8, !tbaa !1398, !dbg !2020
	%153 = getelementptr double, double*  %152, i64  %151, !dbg !2020
	%154 = load double, double*  %153, align 8, !tbaa !1399, !dbg !2020
	store double  %154, double* %..inline.addr.560, align 8, !tbaa !1403, !dbg !2020

	%155 = sub i32  %150, 1, !dbg !2020
	store i32  %155, i32* %..inline.addr.490, align 4, !tbaa !1401, !dbg !2020
	%156 = load i32, i32* %..inline.addr.480, align 4, !tbaa !1401, !dbg !2009
	%157 = sext i32  %156 to i64, !dbg !2009
	%158 = sub i64  %157, 1, !dbg !2009
	store i64  %158, i64* %.TRP0008.addr, align 8, !tbaa !1449, !dbg !2009
	store i64 0, i64* %.ndk0036.addr, align 8, !tbaa !1451, !dbg !2021
	%159 = icmp sle i64  %158, 0, !dbg !2021
	br i1  %159, label %L.B0315, label %L.B0326, !dbg !2021
L.B0326:
	store i64  %158, i64* %.lcr056146.addr, align 8, !tbaa !1451, !dbg !2020
	br label %L.B0314
L.B0314:
	%160 = load i32, i32* %..inline.addr.480, align 4, !tbaa !1401, !dbg !2021
	%161 = load i64, i64* %.ndk0036.addr, align 8, !tbaa !1451, !dbg !2021
	%162 = trunc i64  %161 to i32, !dbg !2021
	%163 = sub i32  %160,  %162, !dbg !2021
	%164 = sext i32  %163 to i64, !dbg !2021
	%165 = sub i64  %164, 1, !dbg !2021
	store i64  %165, i64* %.TRP0009.addr, align 8, !tbaa !1449, !dbg !2021
	store i64 0, i64* %.ndk0038.addr, align 8, !tbaa !1451, !dbg !2021
	%166 = icmp sle i64  %165, 0, !dbg !2021
	br i1  %166, label %L.B0317, label %L.B0327, !dbg !2021
L.B0327:
	store i64  %165, i64* %.lcr056147.addr, align 8, !tbaa !1451, !dbg !2020
	store i32  %162, i32* %.r7.0577.addr, align 4, !tbaa !1401, !dbg !2020
	%167 = bitcast [20 x double]* %..inline.addr.550 to i8*, !dbg !2020
	store i8*  %167, i8** %.G0071.addr, align 8, !tbaa !1398, !dbg !2020
	br label %L.B0316
L.B0316:
	%168 = load i64, i64* %.ndk0038.addr, align 8, !tbaa !1451, !dbg !2021
	%169 = trunc i64  %168 to i32, !dbg !2021
	%170 = sext i32  %169 to i64, !dbg !2021
	%171 = load double*, double** %..inline.addr.460, align 8, !tbaa !1398, !dbg !2021
	%172 = getelementptr double, double*  %171, i64  %170, !dbg !2021
	%173 = load double, double*  %172, align 8, !tbaa !1399, !dbg !2021
	%174 = load i8*, i8** %.G0071.addr, align 8, !tbaa !1398, !dbg !2021
	%175 = getelementptr i8, i8*  %174, i64 8, !dbg !2021
	%176 = bitcast i8*  %175 to double*, !dbg !2021
	%177 = load double, double*  %176, align 8, !tbaa !1399, !dbg !2021
	%178 = bitcast [20 x double]* %..inline.addr.540 to double*, !dbg !2021
	%179 = getelementptr double, double*  %178, i64  %170, !dbg !2021
	%180 = load double, double*  %179, align 8, !tbaa !1399, !dbg !2021
	%181 = fsub double  %177,  %180, !dbg !2021
	%182 = load i32, i32* %.r7.0577.addr, align 4, !tbaa !1401, !dbg !2021
	%183 = add i32  %182,  %169, !dbg !2021
	%184 = sext i32  %183 to i64, !dbg !2021
	%185 = bitcast double*  %171 to i8*, !dbg !2021
	%186 = getelementptr i8, i8*  %185, i64 8, !dbg !2021
	%187 = bitcast i8*  %186 to double*, !dbg !2021
	%188 = getelementptr double, double*  %187, i64  %184, !dbg !2021
	%189 = load double, double*  %188, align 8, !tbaa !1399, !dbg !2021
	%190 = fsub double  %173,  %189, !dbg !2021
	%191 = fdiv double  %181,  %190, !dbg !2021
	%192 = fmul double  %191,  %189, !dbg !2021
	store double  %192, double*  %179, align 8, !tbaa !1399, !dbg !2021
	%193 = fmul double  %173,  %191, !dbg !2021
	%194 = bitcast i8*  %174 to double*, !dbg !2021
	store double  %193, double*  %194, align 8, !tbaa !1399, !dbg !2021
	%195 = add i64  %168, 1, !dbg !2021
	store i64  %195, i64* %.ndk0038.addr, align 8, !tbaa !1451, !dbg !2021
	store i8*  %175, i8** %.G0071.addr, align 8, !tbaa !1398, !dbg !2020
	%196 = load i64, i64* %.lcr056147.addr, align 8, !tbaa !1451, !dbg !2021
	%197 = icmp slt i64  %195,  %196, !dbg !2021
	br i1  %197, label %L.B0316, label %L.B0317, !dbg !2021
L.B0317:
	%198 = load i32, i32* %..inline.addr.490, align 4, !tbaa !1401, !dbg !2021
	%199 = add i32  %198, 1, !dbg !2021
	%200 = mul i32  %199, 2, !dbg !2021
	%201 = load i32, i32* %..inline.addr.480, align 4, !tbaa !1401, !dbg !2021
	%202 = load i64, i64* %.ndk0036.addr, align 8, !tbaa !1451, !dbg !2021
	%203 = trunc i64  %202 to i32, !dbg !2021
	%204 = sub i32  %201,  %203, !dbg !2021
	%205 = sub i32  %204, 1, !dbg !2021
	%206 = icmp sge i32  %200,  %205, !dbg !2021
	br i1  %206, label %L..inline.12246, label %L.B0328, !dbg !2021
L.B0328:
	%207 = sext i32  %198 to i64, !dbg !2021
	%208 = bitcast [20 x double]* %..inline.addr.550 to i8*, !dbg !2021
	%209 = getelementptr i8, i8*  %208, i64 8, !dbg !2021
	%210 = bitcast i8*  %209 to double*, !dbg !2021
	%211 = getelementptr double, double*  %210, i64  %207, !dbg !2021
	%212 = load double, double*  %211, align 8, !tbaa !1399, !dbg !2021
	store double  %212, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2021
	br label %L..inline.12248, !dbg !2021
L..inline.12246:
	%213 = load i32, i32* %..inline.addr.490, align 4, !tbaa !1401, !dbg !2021
	%214 = sext i32  %213 to i64, !dbg !2021
	%215 = bitcast [20 x double]* %..inline.addr.540 to double*, !dbg !2021
	%216 = getelementptr double, double*  %215, i64  %214, !dbg !2021
	%217 = load double, double*  %216, align 8, !tbaa !1399, !dbg !2021
	store double  %217, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2021

	%218 = sub i32  %213, 1, !dbg !2021
	store i32  %218, i32* %..inline.addr.490, align 4, !tbaa !1401, !dbg !2021
	br label %L..inline.12248
L..inline.12248:
	%219 = load double, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2021
	%220 = load double, double* %..inline.addr.560, align 8, !tbaa !1403, !dbg !2021
	%221 = fadd double  %219,  %220, !dbg !2021
	store double  %221, double* %..inline.addr.560, align 8, !tbaa !1403, !dbg !2021
	%222 = load i64, i64* %.ndk0036.addr, align 8, !tbaa !1451, !dbg !2021
	%223 = add i64  %222, 1, !dbg !2021
	store i64  %223, i64* %.ndk0036.addr, align 8, !tbaa !1451, !dbg !2021
	%224 = load i64, i64* %.lcr056146.addr, align 8, !tbaa !1451, !dbg !2021
	%225 = icmp slt i64  %223,  %224, !dbg !2021
	br i1  %225, label %L.B0314, label %L.B0315, !dbg !2021
L.B0315:
	%226 = bitcast [9 x double]* %..inline.addr.270 to i8*, !dbg !2021
	%227 = load i32, i32* %..inline.addr.290, align 4, !tbaa !1401, !dbg !2021
	%228 = load i32, i32* %..inline.addr.450, align 4, !tbaa !1401, !dbg !2021
	%229 = sub i32  %227,  %228, !dbg !2021
	%230 = sext i32  %229 to i64, !dbg !2021
	%231 = mul i64  %230, 8, !dbg !2021
	%232 = getelementptr i8, i8*  %226, i64  %231, !dbg !2021
	%233 = bitcast double** %..inline.addr.680 to i8**, !dbg !2021
	store i8*  %232, i8**  %233, align 8, !tbaa !1398, !dbg !2021
	%234 = bitcast [9 x double]* %..inline.addr.310 to i8*, !dbg !2021
	%235 = mul i64  %230, 8, !dbg !2021
	%236 = getelementptr i8, i8*  %234, i64  %235, !dbg !2021
	%237 = bitcast double** %..inline.addr.690 to i8**, !dbg !2021
	store i8*  %236, i8**  %237, align 8, !tbaa !1398, !dbg !2021
	store i32  %228, i32* %..inline.addr.700, align 4, !tbaa !1401, !dbg !2021
	store i32 0, i32* %..inline.addr.710, align 4, !tbaa !1401, !dbg !2021
	%238 = load double*, double** %..inline.addr.680, align 8, !tbaa !1398, !dbg !2021
	%239 = load double, double*  %238, align 8, !tbaa !1399, !dbg !2021
	%240 = fsub double -0.00000000e+00,  %239, !dbg !2021
	%241 = call double @llvm.fabs.f64 (double  %240), !dbg !2021
	store double  %241, double* %..inline.addr.720, align 8, !tbaa !1403, !dbg !2021
	store i32 0, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	%242 = icmp sle i32  %228, 0, !dbg !2021
	br i1  %242, label %L..inline.12264, label %L.B0329, !dbg !2021
L.B0329:
	store i32  %228, i32* %.lcr016146.addr, align 4, !tbaa !1401, !dbg !2021
	br label %L..inline.12263
L..inline.12263:
	%243 = load i32, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	%244 = sext i32  %243 to i64, !dbg !2021
	%245 = load double*, double** %..inline.addr.680, align 8, !tbaa !1398, !dbg !2021
	%246 = getelementptr double, double*  %245, i64  %244, !dbg !2021
	%247 = load double, double*  %246, align 8, !tbaa !1399, !dbg !2021
	%248 = fsub double -0.00000000e+00,  %247, !dbg !2021
	%249 = call double @llvm.fabs.f64 (double  %248), !dbg !2021
	store double  %249, double* %..inline.addr.740, align 8, !tbaa !1403, !dbg !2021
	%250 = fcmp une double  %249,  0.00000000000000000E+0, !dbg !2021
	br i1  %250, label %L..inline.12266, label %L.B0330, !dbg !2021
L.B0330:
	%251 = load double*, double** %..inline.addr.690, align 8, !tbaa !1398, !dbg !2021
	%252 = getelementptr double, double*  %251, i64  %244, !dbg !2021
	%253 = load double, double*  %252, align 8, !tbaa !1399, !dbg !2021
	store double  %253, double* %..inline.addr.750, align 8, !tbaa !1403, !dbg !2021
	store double  0.00000000000000000E+0, double* %..inline.addr.760, align 8, !tbaa !1403, !dbg !2021
	br label %L..inline.12269, !dbg !2021
L..inline.12266:
	%254 = load double, double* %..inline.addr.740, align 8, !tbaa !1403, !dbg !2021
	%255 = load double, double* %..inline.addr.720, align 8, !tbaa !1403, !dbg !2021
	%256 = fcmp uge double  %254,  %255, !dbg !2021
	br i1  %256, label %L..inline.12271, label %L.B0331, !dbg !2021
L.B0331:
	%257 = load i32, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	store i32  %257, i32* %..inline.addr.710, align 4, !tbaa !1401, !dbg !2021
	store double  %254, double* %..inline.addr.720, align 8, !tbaa !1403, !dbg !2021
	br label %L..inline.12271
L..inline.12271:
	%258 = load i32, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	%259 = sext i32  %258 to i64, !dbg !2021
	%260 = load double*, double** %..inline.addr.690, align 8, !tbaa !1398, !dbg !2021
	%261 = getelementptr double, double*  %260, i64  %259, !dbg !2021
	%262 = load double, double*  %261, align 8, !tbaa !1399, !dbg !2021
	%263 = bitcast [20 x double]* %..inline.addr.770 to double*, !dbg !2021
	%264 = getelementptr double, double*  %263, i64  %259, !dbg !2021
	store double  %262, double*  %264, align 8, !tbaa !1399, !dbg !2021
	%265 = load double, double*  %261, align 8, !tbaa !1399, !dbg !2021
	%266 = fadd double  %265,  1.00000000000000008E-30, !dbg !2021
	%267 = bitcast [20 x double]* %..inline.addr.780 to double*, !dbg !2021
	%268 = getelementptr double, double*  %267, i64  %259, !dbg !2021
	store double  %266, double*  %268, align 8, !tbaa !1399, !dbg !2021

	%269 = add i32  %258, 1, !dbg !2021
	store i32  %269, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	%270 = load i32, i32* %.lcr016146.addr, align 4, !tbaa !1401, !dbg !2021
	%271 = icmp slt i32  %269,  %270, !dbg !2021
	br i1  %271, label %L..inline.12263, label %L..inline.12264, !dbg !2021
L..inline.12264:
	%272 = load i32, i32* %..inline.addr.710, align 4, !tbaa !1401, !dbg !2021

	%273 = sub i32  %272, 1, !dbg !2021
	store i32  %273, i32* %..inline.addr.710, align 4, !tbaa !1401, !dbg !2021
	%274 = sext i32  %272 to i64, !dbg !2021
	%275 = load double*, double** %..inline.addr.690, align 8, !tbaa !1398, !dbg !2021
	%276 = getelementptr double, double*  %275, i64  %274, !dbg !2021
	%277 = load double, double*  %276, align 8, !tbaa !1399, !dbg !2021
	store double  %277, double* %..inline.addr.750, align 8, !tbaa !1403, !dbg !2021
	store i32 1, i32* %..inline.addr.800, align 4, !tbaa !1401, !dbg !2021
	%278 = load i32, i32* %..inline.addr.700, align 4, !tbaa !1401, !dbg !2021
	%279 = icmp sge i32 1,  %278, !dbg !2021
	br i1  %279, label %L..inline.12276, label %L.B0332, !dbg !2021
L.B0332:
	%280 = sub i32  %278, 1, !dbg !2021
	store i32  %280, i32* %.x6146.addr, align 4, !tbaa !1401, !dbg !2021
	br label %L..inline.12275
L..inline.12275:
	store i32 0, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	%281 = load i32, i32* %..inline.addr.700, align 4, !tbaa !1401, !dbg !2021
	%282 = load i32, i32* %..inline.addr.800, align 4, !tbaa !1401, !dbg !2021
	%283 = sub i32  %281,  %282, !dbg !2021
	%284 = icmp sle i32  %283, 0, !dbg !2021
	br i1  %284, label %L..inline.12278, label %L.B0333, !dbg !2021
L.B0333:
	store i32  %283, i32* %.lcr016147.addr, align 4, !tbaa !1401, !dbg !2021
	%285 = bitcast [20 x double]* %..inline.addr.770 to i8*, !dbg !2021
	%286 = getelementptr i8, i8*  %285, i64 8, !dbg !2021
	store i8*  %286, i8** %.G0071.addr, align 8, !tbaa !1398, !dbg !2021
	br label %L..inline.12277
L..inline.12277:
	%287 = load i8*, i8** %.G0071.addr, align 8, !tbaa !1398, !dbg !2021
	%288 = bitcast i8*  %287 to double*, !dbg !2021
	%289 = load double, double*  %288, align 8, !tbaa !1399, !dbg !2021
	%290 = load i32, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	%291 = sext i32  %290 to i64, !dbg !2021
	%292 = bitcast [20 x double]* %..inline.addr.780 to double*, !dbg !2021
	%293 = getelementptr double, double*  %292, i64  %291, !dbg !2021
	%294 = load double, double*  %293, align 8, !tbaa !1399, !dbg !2021
	%295 = fsub double  %289,  %294, !dbg !2021
	store double  %295, double* %..inline.addr.820, align 8, !tbaa !1403, !dbg !2021
	%296 = load double*, double** %..inline.addr.680, align 8, !tbaa !1398, !dbg !2021
	%297 = getelementptr double, double*  %296, i64  %291, !dbg !2021
	%298 = load double, double*  %297, align 8, !tbaa !1399, !dbg !2021
	%299 = fmul double  %294,  %298, !dbg !2021
	%300 = load i32, i32* %..inline.addr.800, align 4, !tbaa !1401, !dbg !2021
	%301 = add i32  %290,  %300, !dbg !2021
	%302 = sext i32  %301 to i64, !dbg !2021
	%303 = getelementptr double, double*  %296, i64  %302, !dbg !2021
	%304 = load double, double*  %303, align 8, !tbaa !1399, !dbg !2021
	%305 = fdiv double  %299,  %304, !dbg !2021
	store double  %305, double* %..inline.addr.830, align 8, !tbaa !1403, !dbg !2021
	%306 = fsub double  %305,  %289, !dbg !2021
	store double  %306, double* %..inline.addr.840, align 8, !tbaa !1403, !dbg !2021
	%307 = fcmp une double  %306,  0.00000000000000000E+0, !dbg !2021
	br i1  %307, label %L..inline.12282, label %L.B0334, !dbg !2021
L.B0334:
	%308 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !2021
	%309 = icmp ne i8*  %308,  null, !dbg !2021
	br i1  %309, label %L..inline.12287, label %L.B0335, !dbg !2021
L.B0335:
	%310 = bitcast %struct._ZSo* @_ZSt4cerr to i8*, !dbg !2021
	%311 = bitcast %struct._ZSo* @_ZSt4cerr to i8**, !dbg !2021
	%312 = load i8*, i8**  %311, align 8, !tbaa !1398, !dbg !2021
	%313 = getelementptr i8, i8*  %312, i64 18446744073709551592, !dbg !2021
	%314 = bitcast i8*  %313 to i64*, !dbg !2021
	%315 = load i64, i64*  %314, align 8, !tbaa !1399, !dbg !2021
	%316 = getelementptr i8, i8*  %310, i64  %315, !dbg !2021
	%317 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.850 to i8**, !dbg !2021
	store i8*  %316, i8**  %317, align 8, !tbaa !1398, !dbg !2021
	%318 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.850, align 8, !tbaa !1398, !dbg !2021
	%319 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %318 to i8*, !dbg !2021
	%320 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.860 to i8**, !dbg !2021
	store i8*  %319, i8**  %320, align 8, !tbaa !1398, !dbg !2021
	%321 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.860, align 8, !tbaa !1398, !dbg !2021
	%322 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %321 to i8*, !dbg !2021
	%323 = getelementptr i8, i8*  %322, i64 32, !dbg !2021
	%324 = bitcast i8*  %323 to i32*, !dbg !2021
	%325 = load i32, i32*  %324, align 4, !tbaa !1399, !dbg !2021
	%326 = or i32  %325, 1, !dbg !2021
	invoke void  @_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate (%struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %318, i32  %326)
		to label %L.B0336
		unwind label %L_T41505256_7635, !dbg !2021
L.B0336:
	br label %L..inline.12308, !dbg !2021
L..inline.12287:
	%327 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !2021
	%328 = call i64  @strlen (i8*  %327) nounwind, !dbg !2021
	%329 = invoke %struct._ZSo*  @_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l (%struct._ZSo* @_ZSt4cerr, i8*  %327, i64  %328)
		to label %L.B0337
		unwind label %L_T41505256_7635, !dbg !2021
L.B0337:
	%330 = bitcast %struct._ZSo*  %329 to i8*, !dbg !2021
	%331 = bitcast %struct._ZSo** %.Q0018.addr to i8**, !dbg !2021
	store i8*  %330, i8**  %331, align 8, !tbaa !1398, !dbg !2021
	br label %L..inline.12308
L..inline.12308:
	call void  @exit (i32 111) nounwind noreturn, !dbg !2021
	br label %L..inline.12282
L..inline.12282:
	%332 = load double, double* %..inline.addr.820, align 8, !tbaa !1403, !dbg !2021
	%333 = load double, double* %..inline.addr.840, align 8, !tbaa !1403, !dbg !2021
	%334 = fdiv double  %332,  %333, !dbg !2021
	%335 = load i8*, i8** %.G0071.addr, align 8, !tbaa !1398, !dbg !2021
	%336 = bitcast i8*  %335 to double*, !dbg !2021
	%337 = load double, double*  %336, align 8, !tbaa !1399, !dbg !2021
	%338 = fmul double  %334,  %337, !dbg !2021
	%339 = load i32, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	%340 = sext i32  %339 to i64, !dbg !2021
	%341 = bitcast [20 x double]* %..inline.addr.780 to double*, !dbg !2021
	%342 = getelementptr double, double*  %341, i64  %340, !dbg !2021
	store double  %338, double*  %342, align 8, !tbaa !1399, !dbg !2021
	%343 = load double, double* %..inline.addr.830, align 8, !tbaa !1403, !dbg !2021
	%344 = fmul double  %343,  %334, !dbg !2021
	%345 = getelementptr i8, i8*  %335, i64 18446744073709551608, !dbg !2021
	%346 = bitcast i8*  %345 to double*, !dbg !2021
	store double  %344, double*  %346, align 8, !tbaa !1399, !dbg !2021

	%347 = add i32  %339, 1, !dbg !2021
	store i32  %347, i32* %..inline.addr.730, align 4, !tbaa !1401, !dbg !2021
	%348 = getelementptr i8, i8*  %335, i64 8, !dbg !2021
	store i8*  %348, i8** %.G0071.addr, align 8, !tbaa !1398, !dbg !2021
	%349 = load i32, i32* %.lcr016147.addr, align 4, !tbaa !1401, !dbg !2021
	%350 = icmp slt i32  %347,  %349, !dbg !2021
	br i1  %350, label %L..inline.12277, label %L..inline.12278, !dbg !2021
L..inline.12278:
	%351 = load i32, i32* %..inline.addr.710, align 4, !tbaa !1401, !dbg !2021
	%352 = add i32  %351, 1, !dbg !2021
	%353 = mul i32  %352, 2, !dbg !2021
	%354 = load i32, i32* %..inline.addr.700, align 4, !tbaa !1401, !dbg !2021
	%355 = load i32, i32* %..inline.addr.800, align 4, !tbaa !1401, !dbg !2021
	%356 = sub i32  %354,  %355, !dbg !2021
	%357 = icmp sge i32  %353,  %356, !dbg !2021
	br i1  %357, label %L..inline.12317, label %L.B0338, !dbg !2021
L.B0338:
	%358 = sext i32  %351 to i64, !dbg !2021
	%359 = bitcast [20 x double]* %..inline.addr.770 to i8*, !dbg !2021
	%360 = getelementptr i8, i8*  %359, i64 8, !dbg !2021
	%361 = bitcast i8*  %360 to double*, !dbg !2021
	%362 = getelementptr double, double*  %361, i64  %358, !dbg !2021
	%363 = load double, double*  %362, align 8, !tbaa !1399, !dbg !2021
	store double  %363, double* %..inline.addr.900, align 8, !tbaa !1403, !dbg !2021
	br label %L..inline.12319, !dbg !2021
L..inline.12317:
	%364 = load i32, i32* %..inline.addr.710, align 4, !tbaa !1401, !dbg !2021

	%365 = sub i32  %364, 1, !dbg !2021
	store i32  %365, i32* %..inline.addr.710, align 4, !tbaa !1401, !dbg !2021
	%366 = sext i32  %364 to i64, !dbg !2021
	%367 = bitcast [20 x double]* %..inline.addr.780 to double*, !dbg !2021
	%368 = getelementptr double, double*  %367, i64  %366, !dbg !2021
	%369 = load double, double*  %368, align 8, !tbaa !1399, !dbg !2021
	store double  %369, double* %..inline.addr.900, align 8, !tbaa !1403, !dbg !2021
	br label %L..inline.12319
L..inline.12319:
	%370 = load double, double* %..inline.addr.900, align 8, !tbaa !1403, !dbg !2021
	store double  %370, double* %..inline.addr.760, align 8, !tbaa !1403, !dbg !2021
	%371 = load double, double* %..inline.addr.750, align 8, !tbaa !1403, !dbg !2021
	%372 = fadd double  %370,  %371, !dbg !2021
	store double  %372, double* %..inline.addr.750, align 8, !tbaa !1403, !dbg !2021
	%373 = load i32, i32* %..inline.addr.800, align 4, !tbaa !1401, !dbg !2021

	%374 = add i32  %373, 1, !dbg !2021
	store i32  %374, i32* %..inline.addr.800, align 4, !tbaa !1401, !dbg !2021
	%375 = load i32, i32* %.x6146.addr, align 4, !tbaa !1401, !dbg !2021
	%376 = sub i32  %375, 1, !dbg !2021
	store i32  %376, i32* %.x6146.addr, align 4, !tbaa !1401, !dbg !2021
	%377 = icmp sgt i32  %376, 0, !dbg !2021
	br i1  %377, label %L..inline.12275, label %L..inline.12276, !dbg !2021
L..inline.12276:
	br label %L..inline.12269
L..inline.12269:
	%378 = load double, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2021
	%379 = call double @llvm.fabs.f64 (double  %378), !dbg !2021
	store double  %379, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2021
	%380 = load double, double* %..inline.addr.760, align 8, !tbaa !1403, !dbg !2022
	%381 = call double @llvm.fabs.f64 (double  %380), !dbg !2022
	store double  %381, double* %..inline.addr.760, align 8, !tbaa !1403, !dbg !2022
	%382 = call double @llvm.maxnum.f64 (double  %379, double  %381), !dbg !2023
	%383 = fcmp uno double  %381,  %381, !dbg !2023
	%384 = select i1  %383, double  %379, double  %382, !dbg !2023
	%385 = load double, double* %..inline.addr.560, align 8, !tbaa !1403, !dbg !2024
	%386 = load double, double* %..inline.addr.750, align 8, !tbaa !1403, !dbg !2024
	%387 = fsub double  %385,  %386, !dbg !2024
	%388 = call double @llvm.fabs.f64 (double  %387), !dbg !2024
	%389 = call double @llvm.maxnum.f64 (double  %384, double  %388), !dbg !2025
	%390 = fcmp uno double  %388,  %388, !dbg !2025
	%391 = select i1  %390, double  %384, double  %389, !dbg !2025
	%392 = load double, double* %..inline.addr.8, align 8, !tbaa !1403, !dbg !2026
	%393 = fcmp uge double  %391,  %392, !dbg !2026
	br i1  %393, label %L..inline.12327, label %L.B0339, !dbg !2026
L.B0339:
	%394 = fadd double  %385,  %386, !dbg !2026
	%395 = fmul double  %394,  5.00000000000000000E-1, !dbg !2026
	store double  %395, double* %..inline.addr.440, align 8, !tbaa !1403, !dbg !2026
	br label %L..inline.12328, !dbg !2027
L..inline.12327:
	br label %L..inline.12205
L..inline.12205:
	%396 = load i8*, i8** %.G0086.addr, align 8, !tbaa !1398, !dbg !2028
	%397 = getelementptr i8, i8*  %396, i64 18446744073709551608, !dbg !2028
	%398 = bitcast i8*  %397 to double*, !dbg !2028
	%399 = load double, double*  %398, align 8, !tbaa !1399, !dbg !2028
	%400 = bitcast i8*  %396 to double*, !dbg !2028
	store double  %399, double*  %400, align 8, !tbaa !1399, !dbg !2028
	%401 = load i8*, i8** %.G0084.addr, align 8, !tbaa !1398, !dbg !2028
	%402 = getelementptr i8, i8*  %401, i64 18446744073709551608, !dbg !2028
	%403 = bitcast i8*  %402 to double*, !dbg !2028
	%404 = load double, double*  %403, align 8, !tbaa !1399, !dbg !2028
	%405 = fmul double  %404,  2.50000000000000000E-1, !dbg !2028
	%406 = bitcast i8*  %401 to double*, !dbg !2028
	store double  %405, double*  %406, align 8, !tbaa !1399, !dbg !2028
	%407 = load i32, i32* %..inline.addr.290, align 4, !tbaa !1401, !dbg !2029

	%408 = add i32  %407, 1, !dbg !2030
	store i32  %408, i32* %..inline.addr.290, align 4, !tbaa !1401, !dbg !2030
	%409 = getelementptr i8, i8*  %396, i64 8, !dbg !2016
	store i8*  %409, i8** %.G0086.addr, align 8, !tbaa !1398, !dbg !2016
	%410 = getelementptr i8, i8*  %401, i64 8, !dbg !2016
	store i8*  %410, i8** %.G0084.addr, align 8, !tbaa !1398, !dbg !2016
	%411 = load i32, i32* %.x6145.addr, align 4, !tbaa !1401, !dbg !2016
	%412 = sub i32  %411, 1, !dbg !2016
	store i32  %412, i32* %.x6145.addr, align 4, !tbaa !1401, !dbg !2016
	%413 = icmp sgt i32  %412, 0, !dbg !2030
	br i1  %413, label %L..inline.12178, label %L.N0023, !dbg !2030
L.N0023:
	br label %L..inline.12328
L..inline.12328:
	%414 = load double, double* %..inline.addr.440, align 8, !tbaa !1403, !dbg !2031
	%415 = bitcast [5 x i32 (...)* (...)*]* @_ZTV10Tintxy_f3D to i8*, !dbg !2009
	%416 = getelementptr i8, i8*  %415, i64 16, !dbg !2009
	%417 = bitcast %struct.Tintxy_f3D* %_T40345952_7635.addr to i8**, !dbg !2009
	store i8*  %416, i8**  %417, align 8, !tbaa !1398, !dbg !2009
	%418 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2009
	%419 = getelementptr i8, i8*  %418, i64 16, !dbg !2009
	store i8*  %419, i8**  %417, align 8, !tbaa !1398, !dbg !2009
	ret double  %414, !dbg !2009
L.B0340:
	br label %L.R0006, !dbg !2009
L_T41505256_7635:
	%420 = landingpad %astruct.dt64
	cleanup
	%421 = extractvalue %astruct.dt64  %420, 0, !dbg !2009
	store i8*  %421, i8** %__caught_object_address.addr, align 1, !tbaa !1398, !dbg !2009
	%422 = extractvalue %astruct.dt64  %420, 1, !dbg !2009
	store i32  %422, i32* %__catch_clause_number.addr, align 1, !tbaa !1399, !dbg !2009
	%423 = bitcast [5 x i32 (...)* (...)*]* @_ZTV10Tintxy_f3D to i8*, !dbg !2009
	%424 = getelementptr i8, i8*  %423, i64 16, !dbg !2009
	%425 = bitcast %struct.Tintxy_f3D* %_T40345952_7635.addr to i8**, !dbg !2009
	store i8*  %424, i8**  %425, align 8, !tbaa !1398, !dbg !2009
	%426 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2009
	%427 = getelementptr i8, i8*  %426, i64 16, !dbg !2009
	store i8*  %427, i8**  %425, align 8, !tbaa !1398, !dbg !2009
	%428 = load i32, i32* %__catch_clause_number.addr, align 4, !tbaa !1399, !dbg !2009
	%429 = load i8*, i8** %__caught_object_address.addr, align 8, !tbaa !1398, !dbg !2009
	%430 = insertvalue %astruct.dt64 undef, i8*  %429, 0, !dbg !2009
	%431 = insertvalue %astruct.dt64  %430, i32  %428, 1, !dbg !2009
	resume %astruct.dt64  %431 , !dbg !2009
L.R0006:
	ret double 0.0
}
define linkonce_odr void @_ZN11T1DFunctionD1Ev(%struct.T1DFunction* %_T40342992_7637.arg) #0 inlinehint !dbg !2039 {
L.entry:
	%_T40342992_7637.addr = alloca %struct.T1DFunction*, align 8

	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %_T40342992_7637.addr, metadata !2043, metadata !1366), !dbg !2040
	store %struct.T1DFunction* %_T40342992_7637.arg, %struct.T1DFunction** %_T40342992_7637.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %_T40342992_7637.addr, metadata !2044, metadata !1366), !dbg !2040
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2045
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2045
	%2 = load %struct.T1DFunction*, %struct.T1DFunction** %_T40342992_7637.addr, align 8, !tbaa !1398, !dbg !2045
	%3 = bitcast %struct.T1DFunction*  %2 to i8**, !dbg !2045
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2045
	ret void, !dbg !2045
}
define linkonce_odr void @_ZN11T1DFunctionD0Ev(%struct.T1DFunction* %_T40342992_7638.arg) #0 inlinehint !dbg !2047 {
L.entry:
	%_T40342992_7638.addr = alloca %struct.T1DFunction*, align 8

	store %struct.T1DFunction* %_T40342992_7638.arg, %struct.T1DFunction** %_T40342992_7638.addr, align 8, !tbaa !1398
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2055
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2055
	%2 = load %struct.T1DFunction*, %struct.T1DFunction** %_T40342992_7638.addr, align 8, !tbaa !1398, !dbg !2055
	%3 = bitcast %struct.T1DFunction*  %2 to i8**, !dbg !2055
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2055
	%4 = bitcast %struct.T1DFunction*  %2 to i8*, !dbg !2055
	call void  @_ZdlPvm (i8*  %4, i64 8) nounwind, !dbg !2055
	ret void, !dbg !2055
}
define linkonce_odr void @_ZN11T1DFunctionD2Ev(%struct.T1DFunction* %_T40342992_7639.arg) #0 inlinehint !dbg !2057 {
L.entry:
	%_T40342992_7639.addr = alloca %struct.T1DFunction*, align 8

	store %struct.T1DFunction* %_T40342992_7639.arg, %struct.T1DFunction** %_T40342992_7639.addr, align 8, !tbaa !1398
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2065
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2065
	%2 = load %struct.T1DFunction*, %struct.T1DFunction** %_T40342992_7639.addr, align 8, !tbaa !1398, !dbg !2065
	%3 = bitcast %struct.T1DFunction*  %2 to i8**, !dbg !2065
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2065
	ret void, !dbg !2065
}
define linkonce_odr void @_ZN11T1DFunctionC1Ev(%struct.T1DFunction* %_T40342992_7640.arg) #0 inlinehint !dbg !2067 {
L.entry:
	%_T40342992_7640.addr = alloca %struct.T1DFunction*, align 8

	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %_T40342992_7640.addr, metadata !2071, metadata !1366), !dbg !2068
	store %struct.T1DFunction* %_T40342992_7640.arg, %struct.T1DFunction** %_T40342992_7640.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T1DFunction** %_T40342992_7640.addr, metadata !2072, metadata !1366), !dbg !2068
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2073
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2073
	%2 = load %struct.T1DFunction*, %struct.T1DFunction** %_T40342992_7640.addr, align 8, !tbaa !1398, !dbg !2073
	%3 = bitcast %struct.T1DFunction*  %2 to i8**, !dbg !2073
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2073
	ret void, !dbg !2073
}
define linkonce_odr void @_ZN11T1DFunctionC2Ev(%struct.T1DFunction* %_T40342992_7641.arg) #0 inlinehint !dbg !2075 {
L.entry:
	%_T40342992_7641.addr = alloca %struct.T1DFunction*, align 8

	store %struct.T1DFunction* %_T40342992_7641.arg, %struct.T1DFunction** %_T40342992_7641.addr, align 8, !tbaa !1398
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2083
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2083
	%2 = load %struct.T1DFunction*, %struct.T1DFunction** %_T40342992_7641.addr, align 8, !tbaa !1398, !dbg !2083
	%3 = bitcast %struct.T1DFunction*  %2 to i8**, !dbg !2083
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2083
	ret void, !dbg !2083
}
define linkonce_odr void @_ZN11T2DFunctionD1Ev(%struct.T2DFunction* %_T40342992_7642.arg) #0 inlinehint !dbg !2088 {
L.entry:
	%_T40342992_7642.addr = alloca %struct.T2DFunction*, align 8

	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %_T40342992_7642.addr, metadata !2092, metadata !1366), !dbg !2089
	store %struct.T2DFunction* %_T40342992_7642.arg, %struct.T2DFunction** %_T40342992_7642.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %_T40342992_7642.addr, metadata !2093, metadata !1366), !dbg !2089
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2094
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2094
	%2 = load %struct.T2DFunction*, %struct.T2DFunction** %_T40342992_7642.addr, align 8, !tbaa !1398, !dbg !2094
	%3 = bitcast %struct.T2DFunction*  %2 to i8**, !dbg !2094
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2094
	ret void, !dbg !2094
}
define linkonce_odr void @_ZN11T2DFunctionD0Ev(%struct.T2DFunction* %_T40342992_7643.arg) #0 inlinehint !dbg !2098 {
L.entry:
	%_T40342992_7643.addr = alloca %struct.T2DFunction*, align 8

	store %struct.T2DFunction* %_T40342992_7643.arg, %struct.T2DFunction** %_T40342992_7643.addr, align 8, !tbaa !1398
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2106
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2106
	%2 = load %struct.T2DFunction*, %struct.T2DFunction** %_T40342992_7643.addr, align 8, !tbaa !1398, !dbg !2106
	%3 = bitcast %struct.T2DFunction*  %2 to i8**, !dbg !2106
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2106
	%4 = bitcast %struct.T2DFunction*  %2 to i8*, !dbg !2106
	call void  @_ZdlPvm (i8*  %4, i64 8) nounwind, !dbg !2106
	ret void, !dbg !2106
}
define linkonce_odr void @_ZN11T2DFunctionD2Ev(%struct.T2DFunction* %_T40342992_7644.arg) #0 inlinehint !dbg !2108 {
L.entry:
	%_T40342992_7644.addr = alloca %struct.T2DFunction*, align 8

	store %struct.T2DFunction* %_T40342992_7644.arg, %struct.T2DFunction** %_T40342992_7644.addr, align 8, !tbaa !1398
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2116
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2116
	%2 = load %struct.T2DFunction*, %struct.T2DFunction** %_T40342992_7644.addr, align 8, !tbaa !1398, !dbg !2116
	%3 = bitcast %struct.T2DFunction*  %2 to i8**, !dbg !2116
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2116
	ret void, !dbg !2116
}
define linkonce_odr void @_ZN11T2DFunctionC1Ev(%struct.T2DFunction* %_T40342992_7645.arg) #0 inlinehint !dbg !2118 {
L.entry:
	%_T40342992_7645.addr = alloca %struct.T2DFunction*, align 8

	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %_T40342992_7645.addr, metadata !2122, metadata !1366), !dbg !2119
	store %struct.T2DFunction* %_T40342992_7645.arg, %struct.T2DFunction** %_T40342992_7645.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %_T40342992_7645.addr, metadata !2123, metadata !1366), !dbg !2119
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2124
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2124
	%2 = load %struct.T2DFunction*, %struct.T2DFunction** %_T40342992_7645.addr, align 8, !tbaa !1398, !dbg !2124
	%3 = bitcast %struct.T2DFunction*  %2 to i8**, !dbg !2124
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2124
	ret void, !dbg !2124
}
define linkonce_odr void @_ZN11T2DFunctionC2Ev(%struct.T2DFunction* %_T40342992_7646.arg) #0 inlinehint !dbg !2126 {
L.entry:
	%_T40342992_7646.addr = alloca %struct.T2DFunction*, align 8

	store %struct.T2DFunction* %_T40342992_7646.arg, %struct.T2DFunction** %_T40342992_7646.addr, align 8, !tbaa !1398
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2134
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2134
	%2 = load %struct.T2DFunction*, %struct.T2DFunction** %_T40342992_7646.addr, align 8, !tbaa !1398, !dbg !2134
	%3 = bitcast %struct.T2DFunction*  %2 to i8**, !dbg !2134
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2134
	ret void, !dbg !2134
}

%struct.T2D_fix1 = type <{ %struct.T1DFunction, %struct.T2DFunction*, double}> 

define linkonce_odr void @_ZN8T2D_fix1C1ERK11T2DFunctiond(%struct.T2D_fix1* %_T40342992_7647.arg, %struct.T2DFunction* %f1.arg, double %x1.arg) #0 inlinehint !dbg !2138 {
L.entry:
	%_T40342992_7647.addr = alloca %struct.T2D_fix1*, align 8
	%f1.addr = alloca %struct.T2DFunction*, align 8
	%x1.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.T2D_fix1** %_T40342992_7647.addr, metadata !2150, metadata !1366), !dbg !2139
	store %struct.T2D_fix1* %_T40342992_7647.arg, %struct.T2D_fix1** %_T40342992_7647.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T2D_fix1** %_T40342992_7647.addr, metadata !2151, metadata !1366), !dbg !2139
	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %f1.addr, metadata !2152, metadata !1366), !dbg !2139
	store %struct.T2DFunction* %f1.arg, %struct.T2DFunction** %f1.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %f1.addr, metadata !2153, metadata !1366), !dbg !2139
	call void @llvm.dbg.declare (metadata double* %x1.addr, metadata !2154, metadata !1366), !dbg !2139
	store double %x1.arg, double* %x1.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %x1.addr, metadata !2155, metadata !1366), !dbg !2139
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2156
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2156
	%2 = load %struct.T2D_fix1*, %struct.T2D_fix1** %_T40342992_7647.addr, align 8, !tbaa !1398, !dbg !2156
	%3 = bitcast %struct.T2D_fix1*  %2 to i8**, !dbg !2156
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2156
	%4 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T2D_fix1 to i8*, !dbg !2156
	%5 = getelementptr i8, i8*  %4, i64 16, !dbg !2156
	store i8*  %5, i8**  %3, align 8, !tbaa !1398, !dbg !2156
	%6 = load %struct.T2DFunction*, %struct.T2DFunction** %f1.addr, align 8, !tbaa !1398, !dbg !2156
	%7 = bitcast %struct.T2DFunction*  %6 to i8*, !dbg !2156
	%8 = bitcast %struct.T2D_fix1*  %2 to i8*, !dbg !2156
	%9 = getelementptr i8, i8*  %8, i64 8, !dbg !2156
	%10 = bitcast i8*  %9 to i8**, !dbg !2156
	store i8*  %7, i8**  %10, align 8, !tbaa !1398, !dbg !2156
	%11 = load double, double* %x1.addr, align 8, !tbaa !1403, !dbg !2156
	%12 = getelementptr i8, i8*  %8, i64 16, !dbg !2156
	%13 = bitcast i8*  %12 to double*, !dbg !2156
	store double  %11, double*  %13, align 8, !tbaa !1399, !dbg !2156
	ret void, !dbg !2156
}
define linkonce_odr void @_ZN8T2D_fix1C2ERK11T2DFunctiond(%struct.T2D_fix1* %_T40342992_7648.arg, %struct.T2DFunction* %_T40343288_7648.arg, double %_T40343584_7648.arg) #0 inlinehint !dbg !2160 {
L.entry:
	%_T40342992_7648.addr = alloca %struct.T2D_fix1*, align 8
	%_T40343288_7648.addr = alloca %struct.T2DFunction*, align 8
	%_T40343584_7648.addr = alloca double, align 8
	%..inline.addr = alloca %struct.T2D_fix1*, align 8
	%..inline.addr.1 = alloca %struct.T2DFunction*, align 8
	%..inline.addr.2 = alloca double, align 8

	store %struct.T2D_fix1* %_T40342992_7648.arg, %struct.T2D_fix1** %_T40342992_7648.addr, align 8, !tbaa !1398
	store %struct.T2DFunction* %_T40343288_7648.arg, %struct.T2DFunction** %_T40343288_7648.addr, align 8, !tbaa !1398
	store double %_T40343584_7648.arg, double* %_T40343584_7648.addr, align 8, !tbaa !1399
	%0 = load %struct.T2D_fix1*, %struct.T2D_fix1** %_T40342992_7648.addr, align 8, !tbaa !1398, !dbg !2176
	%1 = bitcast %struct.T2D_fix1*  %0 to i8*, !dbg !2176
	%2 = bitcast %struct.T2D_fix1** %..inline.addr to i8**, !dbg !2176
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2176
	%3 = load %struct.T2DFunction*, %struct.T2DFunction** %_T40343288_7648.addr, align 8, !tbaa !1398, !dbg !2176
	%4 = bitcast %struct.T2DFunction*  %3 to i8*, !dbg !2176
	%5 = bitcast %struct.T2DFunction** %..inline.addr.1 to i8**, !dbg !2176
	store i8*  %4, i8**  %5, align 8, !tbaa !1398, !dbg !2176
	%6 = load double, double* %_T40343584_7648.addr, align 8, !tbaa !1403, !dbg !2176
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2176
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2176
	%9 = bitcast %struct.T2D_fix1*  %0 to i8**, !dbg !2176
	store i8*  %8, i8**  %9, align 8, !tbaa !1398, !dbg !2176
	%10 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T2D_fix1 to i8*, !dbg !2176
	%11 = getelementptr i8, i8*  %10, i64 16, !dbg !2176
	%12 = load %struct.T2D_fix1*, %struct.T2D_fix1** %..inline.addr, align 8, !tbaa !1398, !dbg !2176
	%13 = bitcast %struct.T2D_fix1*  %12 to i8**, !dbg !2176
	store i8*  %11, i8**  %13, align 8, !tbaa !1398, !dbg !2176
	%14 = load %struct.T2DFunction*, %struct.T2DFunction** %..inline.addr.1, align 8, !tbaa !1398, !dbg !2176
	%15 = bitcast %struct.T2DFunction*  %14 to i8*, !dbg !2176
	%16 = bitcast %struct.T2D_fix1*  %12 to i8*, !dbg !2176
	%17 = getelementptr i8, i8*  %16, i64 8, !dbg !2176
	%18 = bitcast i8*  %17 to i8**, !dbg !2176
	store i8*  %15, i8**  %18, align 8, !tbaa !1398, !dbg !2176
	%19 = getelementptr i8, i8*  %16, i64 16, !dbg !2176
	%20 = bitcast i8*  %19 to double*, !dbg !2176
	store double  %6, double*  %20, align 8, !tbaa !1399, !dbg !2176
	ret void, !dbg !2176
}
define linkonce_odr double @_ZNK8T2D_fix14callEd(%struct.T2D_fix1* %_T40343176_7649.arg, double %y.arg) #0 inlinehint !dbg !2178 {
L.entry:
	%_T40343176_7649.addr = alloca %struct.T2D_fix1*, align 8
	%y.addr = alloca double, align 8
	%.Q0019.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.T2D_fix1** %_T40343176_7649.addr, metadata !2182, metadata !1366), !dbg !2179
	store %struct.T2D_fix1* %_T40343176_7649.arg, %struct.T2D_fix1** %_T40343176_7649.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T2D_fix1** %_T40343176_7649.addr, metadata !2183, metadata !1366), !dbg !2179
	call void @llvm.dbg.declare (metadata double* %y.addr, metadata !2184, metadata !1366), !dbg !2179
	store double %y.arg, double* %y.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %y.addr, metadata !2185, metadata !1366), !dbg !2179
	%0 = load %struct.T2D_fix1*, %struct.T2D_fix1** %_T40343176_7649.addr, align 8, !tbaa !1398, !dbg !2186
	%1 = bitcast %struct.T2D_fix1*  %0 to i8*, !dbg !2186
	%2 = getelementptr i8, i8*  %1, i64 8, !dbg !2186
	%3 = bitcast i8*  %2 to %struct.T2DFunction**, !dbg !2186
	%4 = load %struct.T2DFunction*, %struct.T2DFunction**  %3, align 8, !tbaa !1398, !dbg !2186
	%5 = getelementptr i8, i8*  %1, i64 16, !dbg !2186
	%6 = bitcast i8*  %5 to double*, !dbg !2186
	%7 = load double, double*  %6, align 8, !tbaa !1399, !dbg !2186
	%8 = load double, double* %y.addr, align 8, !tbaa !1403, !dbg !2186
	%9 = bitcast i8*  %2 to double (%struct.T2DFunction*, double, double)****, !dbg !2186
	%10 = load double (%struct.T2DFunction*, double, double)***, double (%struct.T2DFunction*, double, double)****  %9, align 8, !tbaa !1398, !dbg !2186
	%11 = load double (%struct.T2DFunction*, double, double)**, double (%struct.T2DFunction*, double, double)***  %10, align 8, !tbaa !1398, !dbg !2186
	%12 = load double (%struct.T2DFunction*, double, double)*, double (%struct.T2DFunction*, double, double)**  %11, align 8, !tbaa !1398, !dbg !2186
	%13 = call double  %12 (%struct.T2DFunction*  %4, double  %7, double  %8), !dbg !2186
	ret double  %13, !dbg !2186
}
define linkonce_odr void @_ZN8T2D_fix1D1Ev(%struct.T2D_fix1* %_T40342992_7650.arg) #0 inlinehint !dbg !2190 {
L.entry:
	%_T40342992_7650.addr = alloca %struct.T2D_fix1*, align 8

	call void @llvm.dbg.declare (metadata %struct.T2D_fix1** %_T40342992_7650.addr, metadata !2202, metadata !1366), !dbg !2191
	store %struct.T2D_fix1* %_T40342992_7650.arg, %struct.T2D_fix1** %_T40342992_7650.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T2D_fix1** %_T40342992_7650.addr, metadata !2203, metadata !1366), !dbg !2191
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T2D_fix1 to i8*, !dbg !2204
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2204
	%2 = load %struct.T2D_fix1*, %struct.T2D_fix1** %_T40342992_7650.addr, align 8, !tbaa !1398, !dbg !2204
	%3 = bitcast %struct.T2D_fix1*  %2 to i8**, !dbg !2204
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2204
	%4 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2204
	%5 = getelementptr i8, i8*  %4, i64 16, !dbg !2204
	store i8*  %5, i8**  %3, align 8, !tbaa !1398, !dbg !2204
	ret void, !dbg !2204
}
define linkonce_odr void @_ZN8T2D_fix1D0Ev(%struct.T2D_fix1* %_T40342992_7651.arg) #0 inlinehint !dbg !2206 {
L.entry:
	%_T40342992_7651.addr = alloca %struct.T2D_fix1*, align 8
	%..inline.addr = alloca %struct.T2D_fix1*, align 8

	store %struct.T2D_fix1* %_T40342992_7651.arg, %struct.T2D_fix1** %_T40342992_7651.addr, align 8, !tbaa !1398
	%0 = load %struct.T2D_fix1*, %struct.T2D_fix1** %_T40342992_7651.addr, align 8, !tbaa !1398, !dbg !2222
	%1 = bitcast %struct.T2D_fix1*  %0 to i8*, !dbg !2222
	%2 = bitcast %struct.T2D_fix1** %..inline.addr to i8**, !dbg !2222
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2222
	%3 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T2D_fix1 to i8*, !dbg !2222
	%4 = getelementptr i8, i8*  %3, i64 16, !dbg !2222
	%5 = load %struct.T2D_fix1*, %struct.T2D_fix1** %..inline.addr, align 8, !tbaa !1398, !dbg !2222
	%6 = bitcast %struct.T2D_fix1*  %5 to i8**, !dbg !2222
	store i8*  %4, i8**  %6, align 8, !tbaa !1398, !dbg !2222
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2222
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2222
	store i8*  %8, i8**  %6, align 8, !tbaa !1398, !dbg !2222
	call void  @_ZdlPvm (i8*  %1, i64 24) nounwind, !dbg !2222
	ret void, !dbg !2222
}
define linkonce_odr void @_ZN8T2D_fix1D2Ev(%struct.T2D_fix1* %_T40342992_7652.arg) #0 inlinehint !dbg !2224 {
L.entry:
	%_T40342992_7652.addr = alloca %struct.T2D_fix1*, align 8
	%..inline.addr = alloca %struct.T2D_fix1*, align 8

	store %struct.T2D_fix1* %_T40342992_7652.arg, %struct.T2D_fix1** %_T40342992_7652.addr, align 8, !tbaa !1398
	%0 = load %struct.T2D_fix1*, %struct.T2D_fix1** %_T40342992_7652.addr, align 8, !tbaa !1398, !dbg !2240
	%1 = bitcast %struct.T2D_fix1*  %0 to i8*, !dbg !2240
	%2 = bitcast %struct.T2D_fix1** %..inline.addr to i8**, !dbg !2240
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2240
	%3 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T2D_fix1 to i8*, !dbg !2240
	%4 = getelementptr i8, i8*  %3, i64 16, !dbg !2240
	%5 = load %struct.T2D_fix1*, %struct.T2D_fix1** %..inline.addr, align 8, !tbaa !1398, !dbg !2240
	%6 = bitcast %struct.T2D_fix1*  %5 to i8**, !dbg !2240
	store i8*  %4, i8**  %6, align 8, !tbaa !1398, !dbg !2240
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2240
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2240
	store i8*  %8, i8**  %6, align 8, !tbaa !1398, !dbg !2240
	ret void, !dbg !2240
}

%struct.T3D_fix3 = type <{ %struct.T2DFunction, %struct.T3DFunction*, double}> 

define linkonce_odr void @_ZN8T3D_fix3C1ERK11T3DFunctiond(%struct.T3D_fix3* %_T40342992_7653.arg, %struct.T3DFunction* %f1.arg, double %z1.arg) #0 inlinehint !dbg !2244 {
L.entry:
	%_T40342992_7653.addr = alloca %struct.T3D_fix3*, align 8
	%f1.addr = alloca %struct.T3DFunction*, align 8
	%z1.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.T3D_fix3** %_T40342992_7653.addr, metadata !2256, metadata !1366), !dbg !2245
	store %struct.T3D_fix3* %_T40342992_7653.arg, %struct.T3D_fix3** %_T40342992_7653.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T3D_fix3** %_T40342992_7653.addr, metadata !2257, metadata !1366), !dbg !2245
	call void @llvm.dbg.declare (metadata %struct.T3DFunction** %f1.addr, metadata !2258, metadata !1366), !dbg !2245
	store %struct.T3DFunction* %f1.arg, %struct.T3DFunction** %f1.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T3DFunction** %f1.addr, metadata !2259, metadata !1366), !dbg !2245
	call void @llvm.dbg.declare (metadata double* %z1.addr, metadata !2260, metadata !1366), !dbg !2245
	store double %z1.arg, double* %z1.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %z1.addr, metadata !2261, metadata !1366), !dbg !2245
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2262
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2262
	%2 = load %struct.T3D_fix3*, %struct.T3D_fix3** %_T40342992_7653.addr, align 8, !tbaa !1398, !dbg !2262
	%3 = bitcast %struct.T3D_fix3*  %2 to i8**, !dbg !2262
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2262
	%4 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T3D_fix3 to i8*, !dbg !2262
	%5 = getelementptr i8, i8*  %4, i64 16, !dbg !2262
	store i8*  %5, i8**  %3, align 8, !tbaa !1398, !dbg !2262
	%6 = load %struct.T3DFunction*, %struct.T3DFunction** %f1.addr, align 8, !tbaa !1398, !dbg !2262
	%7 = bitcast %struct.T3DFunction*  %6 to i8*, !dbg !2262
	%8 = bitcast %struct.T3D_fix3*  %2 to i8*, !dbg !2262
	%9 = getelementptr i8, i8*  %8, i64 8, !dbg !2262
	%10 = bitcast i8*  %9 to i8**, !dbg !2262
	store i8*  %7, i8**  %10, align 8, !tbaa !1398, !dbg !2262
	%11 = load double, double* %z1.addr, align 8, !tbaa !1403, !dbg !2262
	%12 = getelementptr i8, i8*  %8, i64 16, !dbg !2262
	%13 = bitcast i8*  %12 to double*, !dbg !2262
	store double  %11, double*  %13, align 8, !tbaa !1399, !dbg !2262
	ret void, !dbg !2262
}
define linkonce_odr void @_ZN8T3D_fix3C2ERK11T3DFunctiond(%struct.T3D_fix3* %_T40342992_7654.arg, %struct.T3DFunction* %_T40343288_7654.arg, double %_T40343584_7654.arg) #0 inlinehint !dbg !2266 {
L.entry:
	%_T40342992_7654.addr = alloca %struct.T3D_fix3*, align 8
	%_T40343288_7654.addr = alloca %struct.T3DFunction*, align 8
	%_T40343584_7654.addr = alloca double, align 8
	%..inline.addr = alloca %struct.T3D_fix3*, align 8
	%..inline.addr.1 = alloca %struct.T3DFunction*, align 8
	%..inline.addr.2 = alloca double, align 8

	store %struct.T3D_fix3* %_T40342992_7654.arg, %struct.T3D_fix3** %_T40342992_7654.addr, align 8, !tbaa !1398
	store %struct.T3DFunction* %_T40343288_7654.arg, %struct.T3DFunction** %_T40343288_7654.addr, align 8, !tbaa !1398
	store double %_T40343584_7654.arg, double* %_T40343584_7654.addr, align 8, !tbaa !1399
	%0 = load %struct.T3D_fix3*, %struct.T3D_fix3** %_T40342992_7654.addr, align 8, !tbaa !1398, !dbg !2282
	%1 = bitcast %struct.T3D_fix3*  %0 to i8*, !dbg !2282
	%2 = bitcast %struct.T3D_fix3** %..inline.addr to i8**, !dbg !2282
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2282
	%3 = load %struct.T3DFunction*, %struct.T3DFunction** %_T40343288_7654.addr, align 8, !tbaa !1398, !dbg !2282
	%4 = bitcast %struct.T3DFunction*  %3 to i8*, !dbg !2282
	%5 = bitcast %struct.T3DFunction** %..inline.addr.1 to i8**, !dbg !2282
	store i8*  %4, i8**  %5, align 8, !tbaa !1398, !dbg !2282
	%6 = load double, double* %_T40343584_7654.addr, align 8, !tbaa !1403, !dbg !2282
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2282
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2282
	%9 = bitcast %struct.T3D_fix3*  %0 to i8**, !dbg !2282
	store i8*  %8, i8**  %9, align 8, !tbaa !1398, !dbg !2282
	%10 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T3D_fix3 to i8*, !dbg !2282
	%11 = getelementptr i8, i8*  %10, i64 16, !dbg !2282
	%12 = load %struct.T3D_fix3*, %struct.T3D_fix3** %..inline.addr, align 8, !tbaa !1398, !dbg !2282
	%13 = bitcast %struct.T3D_fix3*  %12 to i8**, !dbg !2282
	store i8*  %11, i8**  %13, align 8, !tbaa !1398, !dbg !2282
	%14 = load %struct.T3DFunction*, %struct.T3DFunction** %..inline.addr.1, align 8, !tbaa !1398, !dbg !2282
	%15 = bitcast %struct.T3DFunction*  %14 to i8*, !dbg !2282
	%16 = bitcast %struct.T3D_fix3*  %12 to i8*, !dbg !2282
	%17 = getelementptr i8, i8*  %16, i64 8, !dbg !2282
	%18 = bitcast i8*  %17 to i8**, !dbg !2282
	store i8*  %15, i8**  %18, align 8, !tbaa !1398, !dbg !2282
	%19 = getelementptr i8, i8*  %16, i64 16, !dbg !2282
	%20 = bitcast i8*  %19 to double*, !dbg !2282
	store double  %6, double*  %20, align 8, !tbaa !1399, !dbg !2282
	ret void, !dbg !2282
}
define linkonce_odr double @_ZNK8T3D_fix34callEdd(%struct.T3D_fix3* %_T40343176_7655.arg, double %x.arg, double %y.arg) #0 inlinehint !dbg !2284 {
L.entry:
	%_T40343176_7655.addr = alloca %struct.T3D_fix3*, align 8
	%x.addr = alloca double, align 8
	%y.addr = alloca double, align 8
	%.Q0020.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.T3D_fix3** %_T40343176_7655.addr, metadata !2288, metadata !1366), !dbg !2285
	store %struct.T3D_fix3* %_T40343176_7655.arg, %struct.T3D_fix3** %_T40343176_7655.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T3D_fix3** %_T40343176_7655.addr, metadata !2289, metadata !1366), !dbg !2285
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !2290, metadata !1366), !dbg !2285
	store double %x.arg, double* %x.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !2291, metadata !1366), !dbg !2285
	call void @llvm.dbg.declare (metadata double* %y.addr, metadata !2292, metadata !1366), !dbg !2285
	store double %y.arg, double* %y.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %y.addr, metadata !2293, metadata !1366), !dbg !2285
	%0 = load %struct.T3D_fix3*, %struct.T3D_fix3** %_T40343176_7655.addr, align 8, !tbaa !1398, !dbg !2294
	%1 = bitcast %struct.T3D_fix3*  %0 to i8*, !dbg !2294
	%2 = getelementptr i8, i8*  %1, i64 8, !dbg !2294
	%3 = bitcast i8*  %2 to %struct.T3DFunction**, !dbg !2294
	%4 = load %struct.T3DFunction*, %struct.T3DFunction**  %3, align 8, !tbaa !1398, !dbg !2294
	%5 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !2294
	%6 = load double, double* %y.addr, align 8, !tbaa !1403, !dbg !2294
	%7 = getelementptr i8, i8*  %1, i64 16, !dbg !2294
	%8 = bitcast i8*  %7 to double*, !dbg !2294
	%9 = load double, double*  %8, align 8, !tbaa !1399, !dbg !2294
	%10 = bitcast i8*  %2 to double (%struct.T3DFunction*, double, double, double)****, !dbg !2294
	%11 = load double (%struct.T3DFunction*, double, double, double)***, double (%struct.T3DFunction*, double, double, double)****  %10, align 8, !tbaa !1398, !dbg !2294
	%12 = load double (%struct.T3DFunction*, double, double, double)**, double (%struct.T3DFunction*, double, double, double)***  %11, align 8, !tbaa !1398, !dbg !2294
	%13 = load double (%struct.T3DFunction*, double, double, double)*, double (%struct.T3DFunction*, double, double, double)**  %12, align 8, !tbaa !1398, !dbg !2294
	%14 = call double  %13 (%struct.T3DFunction*  %4, double  %5, double  %6, double  %9), !dbg !2294
	ret double  %14, !dbg !2294
}
define linkonce_odr void @_ZN8T3D_fix3D1Ev(%struct.T3D_fix3* %_T40342992_7656.arg) #0 inlinehint !dbg !2298 {
L.entry:
	%_T40342992_7656.addr = alloca %struct.T3D_fix3*, align 8

	call void @llvm.dbg.declare (metadata %struct.T3D_fix3** %_T40342992_7656.addr, metadata !2310, metadata !1366), !dbg !2299
	store %struct.T3D_fix3* %_T40342992_7656.arg, %struct.T3D_fix3** %_T40342992_7656.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T3D_fix3** %_T40342992_7656.addr, metadata !2311, metadata !1366), !dbg !2299
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T3D_fix3 to i8*, !dbg !2312
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2312
	%2 = load %struct.T3D_fix3*, %struct.T3D_fix3** %_T40342992_7656.addr, align 8, !tbaa !1398, !dbg !2312
	%3 = bitcast %struct.T3D_fix3*  %2 to i8**, !dbg !2312
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2312
	%4 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2312
	%5 = getelementptr i8, i8*  %4, i64 16, !dbg !2312
	store i8*  %5, i8**  %3, align 8, !tbaa !1398, !dbg !2312
	ret void, !dbg !2312
}
define linkonce_odr void @_ZN8T3D_fix3D0Ev(%struct.T3D_fix3* %_T40342992_7657.arg) #0 inlinehint !dbg !2314 {
L.entry:
	%_T40342992_7657.addr = alloca %struct.T3D_fix3*, align 8
	%..inline.addr = alloca %struct.T3D_fix3*, align 8

	store %struct.T3D_fix3* %_T40342992_7657.arg, %struct.T3D_fix3** %_T40342992_7657.addr, align 8, !tbaa !1398
	%0 = load %struct.T3D_fix3*, %struct.T3D_fix3** %_T40342992_7657.addr, align 8, !tbaa !1398, !dbg !2330
	%1 = bitcast %struct.T3D_fix3*  %0 to i8*, !dbg !2330
	%2 = bitcast %struct.T3D_fix3** %..inline.addr to i8**, !dbg !2330
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2330
	%3 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T3D_fix3 to i8*, !dbg !2330
	%4 = getelementptr i8, i8*  %3, i64 16, !dbg !2330
	%5 = load %struct.T3D_fix3*, %struct.T3D_fix3** %..inline.addr, align 8, !tbaa !1398, !dbg !2330
	%6 = bitcast %struct.T3D_fix3*  %5 to i8**, !dbg !2330
	store i8*  %4, i8**  %6, align 8, !tbaa !1398, !dbg !2330
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2330
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2330
	store i8*  %8, i8**  %6, align 8, !tbaa !1398, !dbg !2330
	call void  @_ZdlPvm (i8*  %1, i64 24) nounwind, !dbg !2330
	ret void, !dbg !2330
}
define linkonce_odr void @_ZN8T3D_fix3D2Ev(%struct.T3D_fix3* %_T40342992_7658.arg) #0 inlinehint !dbg !2332 {
L.entry:
	%_T40342992_7658.addr = alloca %struct.T3D_fix3*, align 8
	%..inline.addr = alloca %struct.T3D_fix3*, align 8

	store %struct.T3D_fix3* %_T40342992_7658.arg, %struct.T3D_fix3** %_T40342992_7658.addr, align 8, !tbaa !1398
	%0 = load %struct.T3D_fix3*, %struct.T3D_fix3** %_T40342992_7658.addr, align 8, !tbaa !1398, !dbg !2348
	%1 = bitcast %struct.T3D_fix3*  %0 to i8*, !dbg !2348
	%2 = bitcast %struct.T3D_fix3** %..inline.addr to i8**, !dbg !2348
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2348
	%3 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T3D_fix3 to i8*, !dbg !2348
	%4 = getelementptr i8, i8*  %3, i64 16, !dbg !2348
	%5 = load %struct.T3D_fix3*, %struct.T3D_fix3** %..inline.addr, align 8, !tbaa !1398, !dbg !2348
	%6 = bitcast %struct.T3D_fix3*  %5 to i8**, !dbg !2348
	store i8*  %4, i8**  %6, align 8, !tbaa !1398, !dbg !2348
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2348
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2348
	store i8*  %8, i8**  %6, align 8, !tbaa !1398, !dbg !2348
	ret void, !dbg !2348
}
define linkonce_odr void @_ZN9Tinty_f2DC1ERK11T2DFunctionddd(%struct.Tinty_f2D* %_T41505256_7659.arg, %struct.T2DFunction* %f1.arg, double %ymin.arg, double %ymax.arg, double %absacc.arg) #0 inlinehint !dbg !2352 {
L.entry:
	%_T41505256_7659.addr = alloca %struct.Tinty_f2D*, align 8
	%f1.addr = alloca %struct.T2DFunction*, align 8
	%ymin.addr = alloca double, align 8
	%ymax.addr = alloca double, align 8
	%absacc.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.Tinty_f2D** %_T41505256_7659.addr, metadata !2364, metadata !1366), !dbg !2353
	store %struct.Tinty_f2D* %_T41505256_7659.arg, %struct.Tinty_f2D** %_T41505256_7659.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.Tinty_f2D** %_T41505256_7659.addr, metadata !2365, metadata !1366), !dbg !2353
	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %f1.addr, metadata !2366, metadata !1366), !dbg !2353
	store %struct.T2DFunction* %f1.arg, %struct.T2DFunction** %f1.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T2DFunction** %f1.addr, metadata !2367, metadata !1366), !dbg !2353
	call void @llvm.dbg.declare (metadata double* %ymin.addr, metadata !2368, metadata !1366), !dbg !2353
	store double %ymin.arg, double* %ymin.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %ymin.addr, metadata !2369, metadata !1366), !dbg !2353
	call void @llvm.dbg.declare (metadata double* %ymax.addr, metadata !2370, metadata !1366), !dbg !2353
	store double %ymax.arg, double* %ymax.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %ymax.addr, metadata !2371, metadata !1366), !dbg !2353
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !2372, metadata !1366), !dbg !2353
	store double %absacc.arg, double* %absacc.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !2373, metadata !1366), !dbg !2353
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2374
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2374
	%2 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %_T41505256_7659.addr, align 8, !tbaa !1398, !dbg !2374
	%3 = bitcast %struct.Tinty_f2D*  %2 to i8**, !dbg !2374
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2374
	%4 = bitcast [5 x i32 (...)* (...)*]* @_ZTV9Tinty_f2D to i8*, !dbg !2374
	%5 = getelementptr i8, i8*  %4, i64 16, !dbg !2374
	store i8*  %5, i8**  %3, align 8, !tbaa !1398, !dbg !2374
	%6 = load %struct.T2DFunction*, %struct.T2DFunction** %f1.addr, align 8, !tbaa !1398, !dbg !2374
	%7 = bitcast %struct.T2DFunction*  %6 to i8*, !dbg !2374
	%8 = bitcast %struct.Tinty_f2D*  %2 to i8*, !dbg !2374
	%9 = getelementptr i8, i8*  %8, i64 8, !dbg !2374
	%10 = bitcast i8*  %9 to i8**, !dbg !2374
	store i8*  %7, i8**  %10, align 8, !tbaa !1398, !dbg !2374
	%11 = load double, double* %ymin.addr, align 8, !tbaa !1403, !dbg !2374
	%12 = getelementptr i8, i8*  %8, i64 16, !dbg !2374
	%13 = bitcast i8*  %12 to double*, !dbg !2374
	store double  %11, double*  %13, align 8, !tbaa !1399, !dbg !2374
	%14 = load double, double* %ymax.addr, align 8, !tbaa !1403, !dbg !2374
	%15 = getelementptr i8, i8*  %8, i64 24, !dbg !2374
	%16 = bitcast i8*  %15 to double*, !dbg !2374
	store double  %14, double*  %16, align 8, !tbaa !1399, !dbg !2374
	%17 = load double, double* %absacc.addr, align 8, !tbaa !1403, !dbg !2374
	%18 = getelementptr i8, i8*  %8, i64 32, !dbg !2374
	%19 = bitcast i8*  %18 to double*, !dbg !2374
	store double  %17, double*  %19, align 8, !tbaa !1399, !dbg !2374
	ret void, !dbg !2374
}
define linkonce_odr void @_ZN9Tinty_f2DC2ERK11T2DFunctionddd(%struct.Tinty_f2D* %_T41505256_7660.arg, %struct.T2DFunction* %_T41505552_7660.arg, double %_T41505848_7660.arg, double %_T41506144_7660.arg, double %_T41506440_7660.arg) #0 inlinehint !dbg !2376 {
L.entry:
	%_T41505256_7660.addr = alloca %struct.Tinty_f2D*, align 8
	%_T41505552_7660.addr = alloca %struct.T2DFunction*, align 8
	%_T41505848_7660.addr = alloca double, align 8
	%_T41506144_7660.addr = alloca double, align 8
	%_T41506440_7660.addr = alloca double, align 8
	%..inline.addr = alloca %struct.Tinty_f2D*, align 8
	%..inline.addr.1 = alloca %struct.T2DFunction*, align 8
	%..inline.addr.2 = alloca double, align 8
	%..inline.addr.3 = alloca double, align 8
	%..inline.addr.4 = alloca double, align 8

	store %struct.Tinty_f2D* %_T41505256_7660.arg, %struct.Tinty_f2D** %_T41505256_7660.addr, align 8, !tbaa !1398
	store %struct.T2DFunction* %_T41505552_7660.arg, %struct.T2DFunction** %_T41505552_7660.addr, align 8, !tbaa !1398
	store double %_T41505848_7660.arg, double* %_T41505848_7660.addr, align 8, !tbaa !1399
	store double %_T41506144_7660.arg, double* %_T41506144_7660.addr, align 8, !tbaa !1399
	store double %_T41506440_7660.arg, double* %_T41506440_7660.addr, align 8, !tbaa !1399
	%0 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %_T41505256_7660.addr, align 8, !tbaa !1398, !dbg !2392
	%1 = bitcast %struct.Tinty_f2D*  %0 to i8*, !dbg !2392
	%2 = bitcast %struct.Tinty_f2D** %..inline.addr to i8**, !dbg !2392
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2392
	%3 = load %struct.T2DFunction*, %struct.T2DFunction** %_T41505552_7660.addr, align 8, !tbaa !1398, !dbg !2392
	%4 = bitcast %struct.T2DFunction*  %3 to i8*, !dbg !2392
	%5 = bitcast %struct.T2DFunction** %..inline.addr.1 to i8**, !dbg !2392
	store i8*  %4, i8**  %5, align 8, !tbaa !1398, !dbg !2392
	%6 = load double, double* %_T41505848_7660.addr, align 8, !tbaa !1403, !dbg !2392
	%7 = load double, double* %_T41506144_7660.addr, align 8, !tbaa !1403, !dbg !2392
	%8 = load double, double* %_T41506440_7660.addr, align 8, !tbaa !1403, !dbg !2392
	%9 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2392
	%10 = getelementptr i8, i8*  %9, i64 16, !dbg !2392
	%11 = bitcast %struct.Tinty_f2D*  %0 to i8**, !dbg !2392
	store i8*  %10, i8**  %11, align 8, !tbaa !1398, !dbg !2392
	%12 = bitcast [5 x i32 (...)* (...)*]* @_ZTV9Tinty_f2D to i8*, !dbg !2392
	%13 = getelementptr i8, i8*  %12, i64 16, !dbg !2392
	%14 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %..inline.addr, align 8, !tbaa !1398, !dbg !2392
	%15 = bitcast %struct.Tinty_f2D*  %14 to i8**, !dbg !2392
	store i8*  %13, i8**  %15, align 8, !tbaa !1398, !dbg !2392
	%16 = load %struct.T2DFunction*, %struct.T2DFunction** %..inline.addr.1, align 8, !tbaa !1398, !dbg !2392
	%17 = bitcast %struct.T2DFunction*  %16 to i8*, !dbg !2392
	%18 = bitcast %struct.Tinty_f2D*  %14 to i8*, !dbg !2392
	%19 = getelementptr i8, i8*  %18, i64 8, !dbg !2392
	%20 = bitcast i8*  %19 to i8**, !dbg !2392
	store i8*  %17, i8**  %20, align 8, !tbaa !1398, !dbg !2392
	%21 = getelementptr i8, i8*  %18, i64 16, !dbg !2392
	%22 = bitcast i8*  %21 to double*, !dbg !2392
	store double  %6, double*  %22, align 8, !tbaa !1399, !dbg !2392
	%23 = getelementptr i8, i8*  %18, i64 24, !dbg !2392
	%24 = bitcast i8*  %23 to double*, !dbg !2392
	store double  %7, double*  %24, align 8, !tbaa !1399, !dbg !2392
	%25 = getelementptr i8, i8*  %18, i64 32, !dbg !2392
	%26 = bitcast i8*  %25 to double*, !dbg !2392
	store double  %8, double*  %26, align 8, !tbaa !1399, !dbg !2392
	ret void, !dbg !2392
}
define linkonce_odr double @_ZNK9Tinty_f2D4callEd(%struct.Tinty_f2D* %_T41505256_7661.arg, double %x.arg) #0 inlinehint personality i8* bitcast (i32 (...)* @__gxx_personality_v0 to i8*) !dbg !2394 {
L.entry:
	%_T41505256_7661.addr = alloca %struct.Tinty_f2D*, align 8
	%x.addr = alloca double, align 8
	%..inline.addr = alloca %struct.T2DFunction*, align 8
	%..inline.addr.1 = alloca double, align 8
	%_T41506440_7661.addr = alloca %struct.T2D_fix1, align 8
	%..inline.addr.2 = alloca double, align 8
	%..inline.addr.3 = alloca double, align 8
	%..inline.addr.4 = alloca double, align 8
	%..inline.addr.5 = alloca i32, align 4
	%..inline.addr.6 = alloca [9 x double], align 8
	%..inline.addr.7 = alloca double, align 8
	%..inline.addr.8 = alloca i32, align 4
	%.x31745.addr = alloca i32, align 4
	%..inline.addr.9 = alloca [9 x double], align 8
	%.G0103.addr = alloca i8*, align 8
	%.G0101.addr = alloca i8*, align 8
	%..inline.addr.180 = alloca double, align 8
	%..inline.addr.190 = alloca double, align 8
	%..inline.addr.200 = alloca double*, align 8
	%.Q0021.addr = alloca double, align 8
	%.Q0022.addr = alloca double, align 8
	%..inline.addr.230 = alloca double, align 8
	%..inline.addr.240 = alloca double, align 8
	%..inline.addr.250 = alloca double, align 8
	%.x31746.addr = alloca i32, align 4
	%.Q0023.addr = alloca double, align 8
	%..inline.addr.280 = alloca double, align 8
	%..inline.addr.290 = alloca i32, align 4
	%..inline.addr.300 = alloca double*, align 8
	%..inline.addr.310 = alloca double*, align 8
	%..inline.addr.320 = alloca i32, align 4
	%..inline.addr.330 = alloca i32, align 4
	%..inline.addr.340 = alloca double, align 8
	%.ndi0022.addr = alloca i32, align 4
	%.lcr0131746.addr = alloca i32, align 4
	%..inline.addr.370 = alloca double, align 8
	%..inline.addr.380 = alloca [20 x double], align 8
	%..inline.addr.390 = alloca [20 x double], align 8
	%..inline.addr.400 = alloca double, align 8
	%.TRP0010.addr = alloca i64, align 8
	%.ndk0044.addr = alloca i64, align 8
	%.lcr0531746.addr = alloca i64, align 8
	%.TRP0011.addr = alloca i64, align 8
	%.ndk0046.addr = alloca i64, align 8
	%.lcr0531747.addr = alloca i64, align 8
	%.r32.0553.addr = alloca i32, align 4
	%.G0088.addr = alloca i8*, align 8
	%..inline.addr.490 = alloca double, align 8
	%..inline.addr.500 = alloca double, align 8
	%..inline.addr.510 = alloca double, align 8
	%..inline.addr.520 = alloca double*, align 8
	%..inline.addr.530 = alloca double*, align 8
	%..inline.addr.540 = alloca i32, align 4
	%..inline.addr.550 = alloca i32, align 4
	%..inline.addr.560 = alloca double, align 8
	%..inline.addr.570 = alloca i32, align 4
	%..inline.addr.580 = alloca double, align 8
	%..inline.addr.590 = alloca double, align 8
	%..inline.addr.600 = alloca double, align 8
	%..inline.addr.610 = alloca [20 x double], align 8
	%..inline.addr.620 = alloca [20 x double], align 8
	%.D0046.addr = alloca i32, align 4
	%..inline.addr.640 = alloca i32, align 4
	%.lcr0131747.addr = alloca i32, align 4
	%..inline.addr.660 = alloca double, align 8
	%..inline.addr.670 = alloca double, align 8
	%..inline.addr.680 = alloca double, align 8
	%..inline.addr.690 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.700 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.710 = alloca i32, align 4
	%..inline.addr.720 = alloca i64, align 8
	%.Q0024.addr = alloca %struct._ZSo*, align 8
	%..inline.addr.740 = alloca double, align 8
	%.D0047.addr = alloca i32, align 4
	%..inline.addr.760 = alloca double, align 8
	%..inline.addr.770 = alloca double, align 8
	%..inline.addr.780 = alloca double, align 8
	%_T41506736_7661.addr = alloca double, align 8
	%__caught_object_address.addr = alloca i8*, align 8
	%__catch_clause_number.addr = alloca i32, align 4

	call void @llvm.dbg.declare (metadata %struct.Tinty_f2D** %_T41505256_7661.addr, metadata !2486, metadata !1366), !dbg !2395
	store %struct.Tinty_f2D* %_T41505256_7661.arg, %struct.Tinty_f2D** %_T41505256_7661.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.Tinty_f2D** %_T41505256_7661.addr, metadata !2487, metadata !1366), !dbg !2395
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !2488, metadata !1366), !dbg !2395
	store double %x.arg, double* %x.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %x.addr, metadata !2489, metadata !1366), !dbg !2395
	%0 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %_T41505256_7661.addr, align 8, !tbaa !1398, !dbg !2490
	%1 = bitcast %struct.Tinty_f2D*  %0 to i8*, !dbg !2490
	%2 = getelementptr i8, i8*  %1, i64 8, !dbg !2490
	%3 = bitcast i8*  %2 to i8**, !dbg !2490
	%4 = load i8*, i8**  %3, align 8, !tbaa !1398, !dbg !2490
	%5 = bitcast %struct.T2DFunction** %..inline.addr to i8**, !dbg !2490
	store i8*  %4, i8**  %5, align 8, !tbaa !1398, !dbg !2490
	%6 = load double, double* %x.addr, align 8, !tbaa !1403, !dbg !2490
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2491
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2491
	%9 = bitcast %struct.T2D_fix1* %_T41506440_7661.addr to i8**, !dbg !2491
	store i8*  %8, i8**  %9, align 8, !tbaa !1398, !dbg !2491
	%10 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T2D_fix1 to i8*, !dbg !2491
	%11 = getelementptr i8, i8*  %10, i64 16, !dbg !2491
	store i8*  %11, i8**  %9, align 8, !tbaa !1398, !dbg !2491
	%12 = load %struct.T2DFunction*, %struct.T2DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !2491
	%13 = bitcast %struct.T2DFunction*  %12 to i8*, !dbg !2491
	%14 = bitcast %struct.T2D_fix1* %_T41506440_7661.addr to i8*, !dbg !2491
	%15 = getelementptr i8, i8*  %14, i64 8, !dbg !2491
	%16 = bitcast i8*  %15 to i8**, !dbg !2491
	store i8*  %13, i8**  %16, align 8, !tbaa !1398, !dbg !2491
	%17 = getelementptr i8, i8*  %14, i64 16, !dbg !2491
	%18 = bitcast i8*  %17 to double*, !dbg !2491
	store double  %6, double*  %18, align 8, !tbaa !1399, !dbg !2491
	%19 = getelementptr i8, i8*  %1, i64 16, !dbg !2490
	%20 = bitcast i8*  %19 to double*, !dbg !2490
	%21 = load double, double*  %20, align 8, !tbaa !1399, !dbg !2490
	store double  %21, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !2490
	%22 = getelementptr i8, i8*  %1, i64 24, !dbg !2490
	%23 = bitcast i8*  %22 to double*, !dbg !2490
	%24 = load double, double*  %23, align 8, !tbaa !1399, !dbg !2490
	store double  %24, double* %..inline.addr.3, align 8, !tbaa !1403, !dbg !2490
	%25 = getelementptr i8, i8*  %1, i64 32, !dbg !2490
	%26 = bitcast i8*  %25 to double*, !dbg !2490
	%27 = load double, double*  %26, align 8, !tbaa !1399, !dbg !2490
	store double  %27, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !2490
	store i32 0, i32* %..inline.addr.5, align 4, !tbaa !1401, !dbg !2490
	%28 = bitcast [9 x double]* %..inline.addr.6 to double*, !dbg !2492
	store double  1.00000000000000000E+0, double*  %28, align 8, !tbaa !1399, !dbg !2492
	store double  0.00000000000000000E+0, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2493
	store i32 1, i32* %..inline.addr.8, align 4, !tbaa !1401, !dbg !2494
	store i32 8, i32* %.x31745.addr, align 4, !tbaa !1401, !dbg !2495
	%29 = bitcast [9 x double]* %..inline.addr.9 to i8*, !dbg !2495
	%30 = getelementptr i8, i8*  %29, i64 8, !dbg !2495
	store i8*  %30, i8** %.G0103.addr, align 8, !tbaa !1398, !dbg !2495
	%31 = bitcast [9 x double]* %..inline.addr.6 to i8*, !dbg !2495
	%32 = getelementptr i8, i8*  %31, i64 8, !dbg !2495
	store i8*  %32, i8** %.G0101.addr, align 8, !tbaa !1398, !dbg !2495
	br label %L..inline.13253
L..inline.13253:
	%33 = load double, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !2496
	store double  %33, double* %..inline.addr.180, align 8, !tbaa !1403, !dbg !2496
	%34 = load double, double* %..inline.addr.3, align 8, !tbaa !1403, !dbg !2496
	store double  %34, double* %..inline.addr.190, align 8, !tbaa !1403, !dbg !2496
	%35 = bitcast [9 x double]* %..inline.addr.9 to i8*, !dbg !2496
	%36 = load i32, i32* %..inline.addr.8, align 4, !tbaa !1401, !dbg !2496
	%37 = sext i32  %36 to i64, !dbg !2496
	%38 = sub i64  %37, 1, !dbg !2496
	%39 = mul i64  %38, 8, !dbg !2496
	%40 = getelementptr i8, i8*  %35, i64  %39, !dbg !2496
	%41 = bitcast double** %..inline.addr.200 to i8**, !dbg !2496
	store i8*  %40, i8**  %41, align 8, !tbaa !1398, !dbg !2496
	%42 = icmp ne i32  %36, 1, !dbg !2496
	br i1  %42, label %L..inline.13267, label %L.B0400, !dbg !2496
L.B0400:
	%43 = bitcast %struct.T2D_fix1* %_T41506440_7661.addr to %struct.T1DFunction*, !dbg !2496
	%44 = bitcast %struct.T2D_fix1* %_T41506440_7661.addr to i8**, !dbg !2496
	%45 = load i8*, i8**  %44, align 8, !tbaa !1398, !dbg !2496
	%46 = bitcast i8*  %45 to double (%struct.T1DFunction*, double)**, !dbg !2496
	%47 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %46, align 8, !tbaa !1398, !dbg !2496
	%48 = invoke double  %47 (%struct.T1DFunction*  %43, double  %34)
		to label %L.B0401
		unwind label %L_T41507328_7661, !dbg !2496
L.B0401:
	store double  %48, double* %.Q0021.addr, align 8, !tbaa !1403, !dbg !2496
	%49 = load double, double* %..inline.addr.180, align 8, !tbaa !1403, !dbg !2496
	%50 = load i8*, i8**  %44, align 8, !tbaa !1398, !dbg !2496
	%51 = bitcast i8*  %50 to double (%struct.T1DFunction*, double)**, !dbg !2496
	%52 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %51, align 8, !tbaa !1398, !dbg !2496
	%53 = invoke double  %52 (%struct.T1DFunction*  %43, double  %49)
		to label %L.B0402
		unwind label %L_T41507328_7661, !dbg !2496
L.B0402:
	store double  %53, double* %.Q0022.addr, align 8, !tbaa !1403, !dbg !2496
	%54 = load double, double* %.Q0021.addr, align 8, !tbaa !1403, !dbg !2496
	%55 = fadd double  %54,  %53, !dbg !2496
	%56 = load double, double* %..inline.addr.190, align 8, !tbaa !1403, !dbg !2496
	%57 = load double, double* %..inline.addr.180, align 8, !tbaa !1403, !dbg !2496
	%58 = fsub double  %56,  %57, !dbg !2496
	%59 = fmul double  %58,  5.00000000000000000E-1, !dbg !2496
	%60 = fmul double  %55,  %59, !dbg !2496
	%61 = load double*, double** %..inline.addr.200, align 8, !tbaa !1398, !dbg !2496
	store double  %60, double*  %61, align 8, !tbaa !1399, !dbg !2496
	store i32 1, i32* %..inline.addr.5, align 4, !tbaa !1401, !dbg !2496
	br label %L..inline.13269, !dbg !2496
L..inline.13267:
	%62 = load double, double* %..inline.addr.3, align 8, !tbaa !1403, !dbg !2496
	%63 = load double, double* %..inline.addr.2, align 8, !tbaa !1403, !dbg !2496
	%64 = fsub double  %62,  %63, !dbg !2496
	%65 = load i32, i32* %..inline.addr.5, align 4, !tbaa !1401, !dbg !2496
	%66 = sitofp i32  %65 to double, !dbg !2496
	%67 = fdiv double  %64,  %66, !dbg !2496
	store double  %67, double* %..inline.addr.230, align 8, !tbaa !1403, !dbg !2496
	%68 = call double @llvm.fma.f64 (double  %67, double  5.00000000000000000E-1, double  %63), !dbg !2496
	store double  %68, double* %..inline.addr.240, align 8, !tbaa !1403, !dbg !2496
	store double  0.00000000000000000E+0, double* %..inline.addr.250, align 8, !tbaa !1403, !dbg !2496
	%69 = icmp sle i32  %65, 0, !dbg !2496
	br i1  %69, label %L..inline.13278, label %L.B0403, !dbg !2496
L.B0403:
	store i32  %65, i32* %.x31746.addr, align 4, !tbaa !1401, !dbg !2496
	br label %L..inline.13277
L..inline.13277:
	%70 = bitcast %struct.T2D_fix1* %_T41506440_7661.addr to %struct.T1DFunction*, !dbg !2496
	%71 = load double, double* %..inline.addr.240, align 8, !tbaa !1403, !dbg !2496
	%72 = bitcast %struct.T2D_fix1* %_T41506440_7661.addr to i8**, !dbg !2496
	%73 = load i8*, i8**  %72, align 8, !tbaa !1398, !dbg !2496
	%74 = bitcast i8*  %73 to double (%struct.T1DFunction*, double)**, !dbg !2496
	%75 = load double (%struct.T1DFunction*, double)*, double (%struct.T1DFunction*, double)**  %74, align 8, !tbaa !1398, !dbg !2496
	%76 = invoke double  %75 (%struct.T1DFunction*  %70, double  %71)
		to label %L.B0404
		unwind label %L_T41507328_7661, !dbg !2496
L.B0404:
	store double  %76, double* %.Q0023.addr, align 8, !tbaa !1403, !dbg !2496
	%77 = load double, double* %..inline.addr.250, align 8, !tbaa !1403, !dbg !2496
	%78 = fadd double  %76,  %77, !dbg !2496
	store double  %78, double* %..inline.addr.250, align 8, !tbaa !1403, !dbg !2496
	%79 = load double, double* %..inline.addr.230, align 8, !tbaa !1403, !dbg !2496
	%80 = load double, double* %..inline.addr.240, align 8, !tbaa !1403, !dbg !2496
	%81 = fadd double  %79,  %80, !dbg !2496
	store double  %81, double* %..inline.addr.240, align 8, !tbaa !1403, !dbg !2496
	%82 = load i32, i32* %.x31746.addr, align 4, !tbaa !1401, !dbg !2496
	%83 = sub i32  %82, 1, !dbg !2496
	store i32  %83, i32* %.x31746.addr, align 4, !tbaa !1401, !dbg !2496
	%84 = icmp sgt i32  %83, 0, !dbg !2496
	br i1  %84, label %L..inline.13277, label %L..inline.13278, !dbg !2496
L..inline.13278:
	%85 = load double, double* %..inline.addr.190, align 8, !tbaa !1403, !dbg !2496
	%86 = load double, double* %..inline.addr.180, align 8, !tbaa !1403, !dbg !2496
	%87 = fsub double  %85,  %86, !dbg !2496
	%88 = load double, double* %..inline.addr.250, align 8, !tbaa !1403, !dbg !2496
	%89 = fmul double  %87,  %88, !dbg !2496
	%90 = load i32, i32* %..inline.addr.5, align 4, !tbaa !1401, !dbg !2496
	%91 = sitofp i32  %90 to double, !dbg !2496
	%92 = fdiv double  %89,  %91, !dbg !2496
	%93 = load double*, double** %..inline.addr.200, align 8, !tbaa !1398, !dbg !2496
	%94 = load double, double*  %93, align 8, !tbaa !1399, !dbg !2496
	%95 = fadd double  %92,  %94, !dbg !2496
	%96 = fmul double  %95,  5.00000000000000000E-1, !dbg !2496
	store double  %96, double*  %93, align 8, !tbaa !1399, !dbg !2496
	%97 = mul i32  %90, 2, !dbg !2496
	store i32  %97, i32* %..inline.addr.5, align 4, !tbaa !1401, !dbg !2496
	br label %L..inline.13269
L..inline.13269:
	%98 = load i8*, i8** %.G0103.addr, align 8, !tbaa !1398, !dbg !2496
	%99 = getelementptr i8, i8*  %98, i64 18446744073709551608, !dbg !2496
	%100 = bitcast i8*  %99 to double*, !dbg !2496
	%101 = load double, double*  %100, align 8, !tbaa !1399, !dbg !2496
	store double  %101, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2496
	%102 = load i32, i32* %..inline.addr.8, align 4, !tbaa !1401, !dbg !2497
	%103 = icmp slt i32  %102, 2, !dbg !2497
	br i1  %103, label %L..inline.13280, label %L.B0405, !dbg !2497
L.B0405:
	%104 = icmp slt i32  %102, 5, !dbg !2498
	%105 = select i1  %104, i32  %102, i32 5, !dbg !2498
	store i32  %105, i32* %..inline.addr.290, align 4, !tbaa !1401, !dbg !2498
	%106 = bitcast [9 x double]* %..inline.addr.6 to i8*, !dbg !2499
	%107 = sub i32  %102,  %105, !dbg !2499
	%108 = sext i32  %107 to i64, !dbg !2499
	%109 = mul i64  %108, 8, !dbg !2499
	%110 = getelementptr i8, i8*  %106, i64  %109, !dbg !2499
	%111 = bitcast double** %..inline.addr.300 to i8**, !dbg !2499
	store i8*  %110, i8**  %111, align 8, !tbaa !1398, !dbg !2499
	%112 = bitcast [9 x double]* %..inline.addr.9 to i8*, !dbg !2499
	%113 = mul i64  %108, 8, !dbg !2499
	%114 = getelementptr i8, i8*  %112, i64  %113, !dbg !2499
	%115 = bitcast double** %..inline.addr.310 to i8**, !dbg !2499
	store i8*  %114, i8**  %115, align 8, !tbaa !1398, !dbg !2499
	store i32  %105, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !2499
	store i32 0, i32* %..inline.addr.330, align 4, !tbaa !1401, !dbg !2499
	%116 = load double*, double** %..inline.addr.300, align 8, !tbaa !1398, !dbg !2499
	%117 = load double, double*  %116, align 8, !tbaa !1399, !dbg !2499
	%118 = fsub double -0.00000000e+00,  %117, !dbg !2499
	%119 = call double @llvm.fabs.f64 (double  %118), !dbg !2499
	store double  %119, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !2499
	store i32 0, i32* %.ndi0022.addr, align 4, !tbaa !1401, !dbg !2500
	%120 = icmp slt i32  %102, 5, !dbg !2500
	%121 = select i1  %120, i32  %102, i32 5, !dbg !2500
	%122 = icmp sle i32  %121, 0, !dbg !2500
	br i1  %122, label %L.B0395, label %L.B0406, !dbg !2500
L.B0406:
	%123 = icmp slt i32  %102, 5, !dbg !2499
	%124 = select i1  %123, i32  %102, i32 5, !dbg !2499
	store i32  %124, i32* %.lcr0131746.addr, align 4, !tbaa !1401, !dbg !2499
	br label %L.B0394
L.B0394:
	%125 = load i32, i32* %.ndi0022.addr, align 4, !tbaa !1401, !dbg !2500
	%126 = sext i32  %125 to i64, !dbg !2500
	%127 = load double*, double** %..inline.addr.300, align 8, !tbaa !1398, !dbg !2500
	%128 = getelementptr double, double*  %127, i64  %126, !dbg !2500
	%129 = load double, double*  %128, align 8, !tbaa !1399, !dbg !2500
	%130 = fsub double -0.00000000e+00,  %129, !dbg !2500
	%131 = call double @llvm.fabs.f64 (double  %130), !dbg !2500
	store double  %131, double* %..inline.addr.370, align 8, !tbaa !1403, !dbg !2500
	%132 = load double, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !2500
	%133 = fcmp uge double  %131,  %132, !dbg !2500
	br i1  %133, label %L..inline.13302, label %L.B0407, !dbg !2500
L.B0407:
	store i32  %125, i32* %..inline.addr.330, align 4, !tbaa !1401, !dbg !2500
	store double  %131, double* %..inline.addr.340, align 8, !tbaa !1403, !dbg !2500
	br label %L..inline.13302
L..inline.13302:
	%134 = load i32, i32* %.ndi0022.addr, align 4, !tbaa !1401, !dbg !2500
	%135 = sext i32  %134 to i64, !dbg !2500
	%136 = load double*, double** %..inline.addr.310, align 8, !tbaa !1398, !dbg !2500
	%137 = getelementptr double, double*  %136, i64  %135, !dbg !2500
	%138 = load double, double*  %137, align 8, !tbaa !1399, !dbg !2500
	%139 = bitcast [20 x double]* %..inline.addr.380 to double*, !dbg !2500
	%140 = getelementptr double, double*  %139, i64  %135, !dbg !2500
	store double  %138, double*  %140, align 8, !tbaa !1399, !dbg !2500
	%141 = bitcast [20 x double]* %..inline.addr.390 to double*, !dbg !2500
	%142 = getelementptr double, double*  %141, i64  %135, !dbg !2500
	store double  %138, double*  %142, align 8, !tbaa !1399, !dbg !2500
	%143 = add i32  %134, 1, !dbg !2500
	store i32  %143, i32* %.ndi0022.addr, align 4, !tbaa !1401, !dbg !2500
	%144 = load i32, i32* %.lcr0131746.addr, align 4, !tbaa !1401, !dbg !2500
	%145 = icmp slt i32  %143,  %144, !dbg !2500
	br i1  %145, label %L.B0394, label %L.B0395, !dbg !2500
L.B0395:
	%146 = load i32, i32* %..inline.addr.330, align 4, !tbaa !1401, !dbg !2499
	%147 = sext i32  %146 to i64, !dbg !2499
	%148 = load double*, double** %..inline.addr.310, align 8, !tbaa !1398, !dbg !2499
	%149 = getelementptr double, double*  %148, i64  %147, !dbg !2499
	%150 = load double, double*  %149, align 8, !tbaa !1399, !dbg !2499
	store double  %150, double* %..inline.addr.400, align 8, !tbaa !1403, !dbg !2499

	%151 = sub i32  %146, 1, !dbg !2499
	store i32  %151, i32* %..inline.addr.330, align 4, !tbaa !1401, !dbg !2499
	%152 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !2490
	%153 = sext i32  %152 to i64, !dbg !2490
	%154 = sub i64  %153, 1, !dbg !2490
	store i64  %154, i64* %.TRP0010.addr, align 8, !tbaa !1449, !dbg !2490
	store i64 0, i64* %.ndk0044.addr, align 8, !tbaa !1451, !dbg !2500
	%155 = icmp sle i64  %154, 0, !dbg !2500
	br i1  %155, label %L.B0397, label %L.B0408, !dbg !2500
L.B0408:
	store i64  %154, i64* %.lcr0531746.addr, align 8, !tbaa !1451, !dbg !2499
	br label %L.B0396
L.B0396:
	%156 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !2500
	%157 = load i64, i64* %.ndk0044.addr, align 8, !tbaa !1451, !dbg !2500
	%158 = trunc i64  %157 to i32, !dbg !2500
	%159 = sub i32  %156,  %158, !dbg !2500
	%160 = sext i32  %159 to i64, !dbg !2500
	%161 = sub i64  %160, 1, !dbg !2500
	store i64  %161, i64* %.TRP0011.addr, align 8, !tbaa !1449, !dbg !2500
	store i64 0, i64* %.ndk0046.addr, align 8, !tbaa !1451, !dbg !2500
	%162 = icmp sle i64  %161, 0, !dbg !2500
	br i1  %162, label %L.B0399, label %L.B0409, !dbg !2500
L.B0409:
	store i64  %161, i64* %.lcr0531747.addr, align 8, !tbaa !1451, !dbg !2499
	store i32  %158, i32* %.r32.0553.addr, align 4, !tbaa !1401, !dbg !2499
	%163 = bitcast [20 x double]* %..inline.addr.390 to i8*, !dbg !2499
	store i8*  %163, i8** %.G0088.addr, align 8, !tbaa !1398, !dbg !2499
	br label %L.B0398
L.B0398:
	%164 = load i64, i64* %.ndk0046.addr, align 8, !tbaa !1451, !dbg !2500
	%165 = trunc i64  %164 to i32, !dbg !2500
	%166 = sext i32  %165 to i64, !dbg !2500
	%167 = load double*, double** %..inline.addr.300, align 8, !tbaa !1398, !dbg !2500
	%168 = getelementptr double, double*  %167, i64  %166, !dbg !2500
	%169 = load double, double*  %168, align 8, !tbaa !1399, !dbg !2500
	%170 = load i8*, i8** %.G0088.addr, align 8, !tbaa !1398, !dbg !2500
	%171 = getelementptr i8, i8*  %170, i64 8, !dbg !2500
	%172 = bitcast i8*  %171 to double*, !dbg !2500
	%173 = load double, double*  %172, align 8, !tbaa !1399, !dbg !2500
	%174 = bitcast [20 x double]* %..inline.addr.380 to double*, !dbg !2500
	%175 = getelementptr double, double*  %174, i64  %166, !dbg !2500
	%176 = load double, double*  %175, align 8, !tbaa !1399, !dbg !2500
	%177 = fsub double  %173,  %176, !dbg !2500
	%178 = load i32, i32* %.r32.0553.addr, align 4, !tbaa !1401, !dbg !2500
	%179 = add i32  %178,  %165, !dbg !2500
	%180 = sext i32  %179 to i64, !dbg !2500
	%181 = bitcast double*  %167 to i8*, !dbg !2500
	%182 = getelementptr i8, i8*  %181, i64 8, !dbg !2500
	%183 = bitcast i8*  %182 to double*, !dbg !2500
	%184 = getelementptr double, double*  %183, i64  %180, !dbg !2500
	%185 = load double, double*  %184, align 8, !tbaa !1399, !dbg !2500
	%186 = fsub double  %169,  %185, !dbg !2500
	%187 = fdiv double  %177,  %186, !dbg !2500
	%188 = fmul double  %187,  %185, !dbg !2500
	store double  %188, double*  %175, align 8, !tbaa !1399, !dbg !2500
	%189 = fmul double  %169,  %187, !dbg !2500
	%190 = bitcast i8*  %170 to double*, !dbg !2500
	store double  %189, double*  %190, align 8, !tbaa !1399, !dbg !2500
	%191 = add i64  %164, 1, !dbg !2500
	store i64  %191, i64* %.ndk0046.addr, align 8, !tbaa !1451, !dbg !2500
	store i8*  %171, i8** %.G0088.addr, align 8, !tbaa !1398, !dbg !2499
	%192 = load i64, i64* %.lcr0531747.addr, align 8, !tbaa !1451, !dbg !2500
	%193 = icmp slt i64  %191,  %192, !dbg !2500
	br i1  %193, label %L.B0398, label %L.B0399, !dbg !2500
L.B0399:
	%194 = load i32, i32* %..inline.addr.330, align 4, !tbaa !1401, !dbg !2500
	%195 = add i32  %194, 1, !dbg !2500
	%196 = mul i32  %195, 2, !dbg !2500
	%197 = load i32, i32* %..inline.addr.320, align 4, !tbaa !1401, !dbg !2500
	%198 = load i64, i64* %.ndk0044.addr, align 8, !tbaa !1451, !dbg !2500
	%199 = trunc i64  %198 to i32, !dbg !2500
	%200 = sub i32  %197,  %199, !dbg !2500
	%201 = sub i32  %200, 1, !dbg !2500
	%202 = icmp sge i32  %196,  %201, !dbg !2500
	br i1  %202, label %L..inline.13321, label %L.B0410, !dbg !2500
L.B0410:
	%203 = sext i32  %194 to i64, !dbg !2500
	%204 = bitcast [20 x double]* %..inline.addr.390 to i8*, !dbg !2500
	%205 = getelementptr i8, i8*  %204, i64 8, !dbg !2500
	%206 = bitcast i8*  %205 to double*, !dbg !2500
	%207 = getelementptr double, double*  %206, i64  %203, !dbg !2500
	%208 = load double, double*  %207, align 8, !tbaa !1399, !dbg !2500
	store double  %208, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2500
	br label %L..inline.13323, !dbg !2500
L..inline.13321:
	%209 = load i32, i32* %..inline.addr.330, align 4, !tbaa !1401, !dbg !2500
	%210 = sext i32  %209 to i64, !dbg !2500
	%211 = bitcast [20 x double]* %..inline.addr.380 to double*, !dbg !2500
	%212 = getelementptr double, double*  %211, i64  %210, !dbg !2500
	%213 = load double, double*  %212, align 8, !tbaa !1399, !dbg !2500
	store double  %213, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2500

	%214 = sub i32  %209, 1, !dbg !2500
	store i32  %214, i32* %..inline.addr.330, align 4, !tbaa !1401, !dbg !2500
	br label %L..inline.13323
L..inline.13323:
	%215 = load double, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2500
	%216 = load double, double* %..inline.addr.400, align 8, !tbaa !1403, !dbg !2500
	%217 = fadd double  %215,  %216, !dbg !2500
	store double  %217, double* %..inline.addr.400, align 8, !tbaa !1403, !dbg !2500
	%218 = load i64, i64* %.ndk0044.addr, align 8, !tbaa !1451, !dbg !2500
	%219 = add i64  %218, 1, !dbg !2500
	store i64  %219, i64* %.ndk0044.addr, align 8, !tbaa !1451, !dbg !2500
	%220 = load i64, i64* %.lcr0531746.addr, align 8, !tbaa !1451, !dbg !2500
	%221 = icmp slt i64  %219,  %220, !dbg !2500
	br i1  %221, label %L.B0396, label %L.B0397, !dbg !2500
L.B0397:
	%222 = bitcast [9 x double]* %..inline.addr.6 to i8*, !dbg !2500
	%223 = load i32, i32* %..inline.addr.8, align 4, !tbaa !1401, !dbg !2500
	%224 = load i32, i32* %..inline.addr.290, align 4, !tbaa !1401, !dbg !2500
	%225 = sub i32  %223,  %224, !dbg !2500
	%226 = sext i32  %225 to i64, !dbg !2500
	%227 = mul i64  %226, 8, !dbg !2500
	%228 = getelementptr i8, i8*  %222, i64  %227, !dbg !2500
	%229 = bitcast double** %..inline.addr.520 to i8**, !dbg !2500
	store i8*  %228, i8**  %229, align 8, !tbaa !1398, !dbg !2500
	%230 = bitcast [9 x double]* %..inline.addr.9 to i8*, !dbg !2500
	%231 = mul i64  %226, 8, !dbg !2500
	%232 = getelementptr i8, i8*  %230, i64  %231, !dbg !2500
	%233 = bitcast double** %..inline.addr.530 to i8**, !dbg !2500
	store i8*  %232, i8**  %233, align 8, !tbaa !1398, !dbg !2500
	store i32  %224, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !2500
	store i32 0, i32* %..inline.addr.550, align 4, !tbaa !1401, !dbg !2500
	%234 = load double*, double** %..inline.addr.520, align 8, !tbaa !1398, !dbg !2500
	%235 = load double, double*  %234, align 8, !tbaa !1399, !dbg !2500
	%236 = fsub double -0.00000000e+00,  %235, !dbg !2500
	%237 = call double @llvm.fabs.f64 (double  %236), !dbg !2500
	store double  %237, double* %..inline.addr.560, align 8, !tbaa !1403, !dbg !2500
	store i32 0, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	%238 = icmp sle i32  %224, 0, !dbg !2500
	br i1  %238, label %L..inline.13339, label %L.B0411, !dbg !2500
L.B0411:
	store i32  %224, i32* %.lcr0131746.addr, align 4, !tbaa !1401, !dbg !2500
	br label %L..inline.13338
L..inline.13338:
	%239 = load i32, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	%240 = sext i32  %239 to i64, !dbg !2500
	%241 = load double*, double** %..inline.addr.520, align 8, !tbaa !1398, !dbg !2500
	%242 = getelementptr double, double*  %241, i64  %240, !dbg !2500
	%243 = load double, double*  %242, align 8, !tbaa !1399, !dbg !2500
	%244 = fsub double -0.00000000e+00,  %243, !dbg !2500
	%245 = call double @llvm.fabs.f64 (double  %244), !dbg !2500
	store double  %245, double* %..inline.addr.580, align 8, !tbaa !1403, !dbg !2500
	%246 = fcmp une double  %245,  0.00000000000000000E+0, !dbg !2500
	br i1  %246, label %L..inline.13341, label %L.B0412, !dbg !2500
L.B0412:
	%247 = load double*, double** %..inline.addr.530, align 8, !tbaa !1398, !dbg !2500
	%248 = getelementptr double, double*  %247, i64  %240, !dbg !2500
	%249 = load double, double*  %248, align 8, !tbaa !1399, !dbg !2500
	store double  %249, double* %..inline.addr.590, align 8, !tbaa !1403, !dbg !2500
	store double  0.00000000000000000E+0, double* %..inline.addr.600, align 8, !tbaa !1403, !dbg !2500
	br label %L..inline.13344, !dbg !2500
L..inline.13341:
	%250 = load double, double* %..inline.addr.580, align 8, !tbaa !1403, !dbg !2500
	%251 = load double, double* %..inline.addr.560, align 8, !tbaa !1403, !dbg !2500
	%252 = fcmp uge double  %250,  %251, !dbg !2500
	br i1  %252, label %L..inline.13346, label %L.B0413, !dbg !2500
L.B0413:
	%253 = load i32, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	store i32  %253, i32* %..inline.addr.550, align 4, !tbaa !1401, !dbg !2500
	store double  %250, double* %..inline.addr.560, align 8, !tbaa !1403, !dbg !2500
	br label %L..inline.13346
L..inline.13346:
	%254 = load i32, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	%255 = sext i32  %254 to i64, !dbg !2500
	%256 = load double*, double** %..inline.addr.530, align 8, !tbaa !1398, !dbg !2500
	%257 = getelementptr double, double*  %256, i64  %255, !dbg !2500
	%258 = load double, double*  %257, align 8, !tbaa !1399, !dbg !2500
	%259 = bitcast [20 x double]* %..inline.addr.610 to double*, !dbg !2500
	%260 = getelementptr double, double*  %259, i64  %255, !dbg !2500
	store double  %258, double*  %260, align 8, !tbaa !1399, !dbg !2500
	%261 = load double, double*  %257, align 8, !tbaa !1399, !dbg !2500
	%262 = fadd double  %261,  1.00000000000000008E-30, !dbg !2500
	%263 = bitcast [20 x double]* %..inline.addr.620 to double*, !dbg !2500
	%264 = getelementptr double, double*  %263, i64  %255, !dbg !2500
	store double  %262, double*  %264, align 8, !tbaa !1399, !dbg !2500

	%265 = add i32  %254, 1, !dbg !2500
	store i32  %265, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	%266 = load i32, i32* %.lcr0131746.addr, align 4, !tbaa !1401, !dbg !2500
	%267 = icmp slt i32  %265,  %266, !dbg !2500
	br i1  %267, label %L..inline.13338, label %L..inline.13339, !dbg !2500
L..inline.13339:
	%268 = load i32, i32* %..inline.addr.550, align 4, !tbaa !1401, !dbg !2500

	%269 = sub i32  %268, 1, !dbg !2500
	store i32  %269, i32* %..inline.addr.550, align 4, !tbaa !1401, !dbg !2500
	%270 = sext i32  %268 to i64, !dbg !2500
	%271 = load double*, double** %..inline.addr.530, align 8, !tbaa !1398, !dbg !2500
	%272 = getelementptr double, double*  %271, i64  %270, !dbg !2500
	%273 = load double, double*  %272, align 8, !tbaa !1399, !dbg !2500
	store double  %273, double* %..inline.addr.590, align 8, !tbaa !1403, !dbg !2500
	store i32 1, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !2500
	%274 = load i32, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !2500
	%275 = icmp sge i32 1,  %274, !dbg !2500
	br i1  %275, label %L..inline.13351, label %L.B0414, !dbg !2500
L.B0414:
	%276 = sub i32  %274, 1, !dbg !2500
	store i32  %276, i32* %.x31746.addr, align 4, !tbaa !1401, !dbg !2500
	br label %L..inline.13350
L..inline.13350:
	store i32 0, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	%277 = load i32, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !2500
	%278 = load i32, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !2500
	%279 = sub i32  %277,  %278, !dbg !2500
	%280 = icmp sle i32  %279, 0, !dbg !2500
	br i1  %280, label %L..inline.13353, label %L.B0415, !dbg !2500
L.B0415:
	store i32  %279, i32* %.lcr0131747.addr, align 4, !tbaa !1401, !dbg !2500
	%281 = bitcast [20 x double]* %..inline.addr.610 to i8*, !dbg !2500
	%282 = getelementptr i8, i8*  %281, i64 8, !dbg !2500
	store i8*  %282, i8** %.G0088.addr, align 8, !tbaa !1398, !dbg !2500
	br label %L..inline.13352
L..inline.13352:
	%283 = load i8*, i8** %.G0088.addr, align 8, !tbaa !1398, !dbg !2500
	%284 = bitcast i8*  %283 to double*, !dbg !2500
	%285 = load double, double*  %284, align 8, !tbaa !1399, !dbg !2500
	%286 = load i32, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	%287 = sext i32  %286 to i64, !dbg !2500
	%288 = bitcast [20 x double]* %..inline.addr.620 to double*, !dbg !2500
	%289 = getelementptr double, double*  %288, i64  %287, !dbg !2500
	%290 = load double, double*  %289, align 8, !tbaa !1399, !dbg !2500
	%291 = fsub double  %285,  %290, !dbg !2500
	store double  %291, double* %..inline.addr.660, align 8, !tbaa !1403, !dbg !2500
	%292 = load double*, double** %..inline.addr.520, align 8, !tbaa !1398, !dbg !2500
	%293 = getelementptr double, double*  %292, i64  %287, !dbg !2500
	%294 = load double, double*  %293, align 8, !tbaa !1399, !dbg !2500
	%295 = fmul double  %290,  %294, !dbg !2500
	%296 = load i32, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !2500
	%297 = add i32  %286,  %296, !dbg !2500
	%298 = sext i32  %297 to i64, !dbg !2500
	%299 = getelementptr double, double*  %292, i64  %298, !dbg !2500
	%300 = load double, double*  %299, align 8, !tbaa !1399, !dbg !2500
	%301 = fdiv double  %295,  %300, !dbg !2500
	store double  %301, double* %..inline.addr.670, align 8, !tbaa !1403, !dbg !2500
	%302 = fsub double  %301,  %285, !dbg !2500
	store double  %302, double* %..inline.addr.680, align 8, !tbaa !1403, !dbg !2500
	%303 = fcmp une double  %302,  0.00000000000000000E+0, !dbg !2500
	br i1  %303, label %L..inline.13357, label %L.B0416, !dbg !2500
L.B0416:
	%304 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !2500
	%305 = icmp ne i8*  %304,  null, !dbg !2500
	br i1  %305, label %L..inline.13362, label %L.B0417, !dbg !2500
L.B0417:
	%306 = bitcast %struct._ZSo* @_ZSt4cerr to i8*, !dbg !2500
	%307 = bitcast %struct._ZSo* @_ZSt4cerr to i8**, !dbg !2500
	%308 = load i8*, i8**  %307, align 8, !tbaa !1398, !dbg !2500
	%309 = getelementptr i8, i8*  %308, i64 18446744073709551592, !dbg !2500
	%310 = bitcast i8*  %309 to i64*, !dbg !2500
	%311 = load i64, i64*  %310, align 8, !tbaa !1399, !dbg !2500
	%312 = getelementptr i8, i8*  %306, i64  %311, !dbg !2500
	%313 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.690 to i8**, !dbg !2500
	store i8*  %312, i8**  %313, align 8, !tbaa !1398, !dbg !2500
	%314 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.690, align 8, !tbaa !1398, !dbg !2500
	%315 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %314 to i8*, !dbg !2500
	%316 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.700 to i8**, !dbg !2500
	store i8*  %315, i8**  %316, align 8, !tbaa !1398, !dbg !2500
	%317 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.700, align 8, !tbaa !1398, !dbg !2500
	%318 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %317 to i8*, !dbg !2500
	%319 = getelementptr i8, i8*  %318, i64 32, !dbg !2500
	%320 = bitcast i8*  %319 to i32*, !dbg !2500
	%321 = load i32, i32*  %320, align 4, !tbaa !1399, !dbg !2500
	%322 = or i32  %321, 1, !dbg !2500
	invoke void  @_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate (%struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %314, i32  %322)
		to label %L.B0418
		unwind label %L_T41507328_7661, !dbg !2500
L.B0418:
	br label %L..inline.13383, !dbg !2500
L..inline.13362:
	%323 = bitcast [21 x i8]* @.S08003 to i8*, !dbg !2500
	%324 = call i64  @strlen (i8*  %323) nounwind, !dbg !2500
	%325 = invoke %struct._ZSo*  @_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l (%struct._ZSo* @_ZSt4cerr, i8*  %323, i64  %324)
		to label %L.B0419
		unwind label %L_T41507328_7661, !dbg !2500
L.B0419:
	%326 = bitcast %struct._ZSo*  %325 to i8*, !dbg !2500
	%327 = bitcast %struct._ZSo** %.Q0024.addr to i8**, !dbg !2500
	store i8*  %326, i8**  %327, align 8, !tbaa !1398, !dbg !2500
	br label %L..inline.13383
L..inline.13383:
	call void  @exit (i32 111) nounwind noreturn, !dbg !2500
	br label %L..inline.13357
L..inline.13357:
	%328 = load double, double* %..inline.addr.660, align 8, !tbaa !1403, !dbg !2500
	%329 = load double, double* %..inline.addr.680, align 8, !tbaa !1403, !dbg !2500
	%330 = fdiv double  %328,  %329, !dbg !2500
	%331 = load i8*, i8** %.G0088.addr, align 8, !tbaa !1398, !dbg !2500
	%332 = bitcast i8*  %331 to double*, !dbg !2500
	%333 = load double, double*  %332, align 8, !tbaa !1399, !dbg !2500
	%334 = fmul double  %330,  %333, !dbg !2500
	%335 = load i32, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	%336 = sext i32  %335 to i64, !dbg !2500
	%337 = bitcast [20 x double]* %..inline.addr.620 to double*, !dbg !2500
	%338 = getelementptr double, double*  %337, i64  %336, !dbg !2500
	store double  %334, double*  %338, align 8, !tbaa !1399, !dbg !2500
	%339 = load double, double* %..inline.addr.670, align 8, !tbaa !1403, !dbg !2500
	%340 = fmul double  %339,  %330, !dbg !2500
	%341 = getelementptr i8, i8*  %331, i64 18446744073709551608, !dbg !2500
	%342 = bitcast i8*  %341 to double*, !dbg !2500
	store double  %340, double*  %342, align 8, !tbaa !1399, !dbg !2500

	%343 = add i32  %335, 1, !dbg !2500
	store i32  %343, i32* %..inline.addr.570, align 4, !tbaa !1401, !dbg !2500
	%344 = getelementptr i8, i8*  %331, i64 8, !dbg !2500
	store i8*  %344, i8** %.G0088.addr, align 8, !tbaa !1398, !dbg !2500
	%345 = load i32, i32* %.lcr0131747.addr, align 4, !tbaa !1401, !dbg !2500
	%346 = icmp slt i32  %343,  %345, !dbg !2500
	br i1  %346, label %L..inline.13352, label %L..inline.13353, !dbg !2500
L..inline.13353:
	%347 = load i32, i32* %..inline.addr.550, align 4, !tbaa !1401, !dbg !2500
	%348 = add i32  %347, 1, !dbg !2500
	%349 = mul i32  %348, 2, !dbg !2500
	%350 = load i32, i32* %..inline.addr.540, align 4, !tbaa !1401, !dbg !2500
	%351 = load i32, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !2500
	%352 = sub i32  %350,  %351, !dbg !2500
	%353 = icmp sge i32  %349,  %352, !dbg !2500
	br i1  %353, label %L..inline.13392, label %L.B0420, !dbg !2500
L.B0420:
	%354 = sext i32  %347 to i64, !dbg !2500
	%355 = bitcast [20 x double]* %..inline.addr.610 to i8*, !dbg !2500
	%356 = getelementptr i8, i8*  %355, i64 8, !dbg !2500
	%357 = bitcast i8*  %356 to double*, !dbg !2500
	%358 = getelementptr double, double*  %357, i64  %354, !dbg !2500
	%359 = load double, double*  %358, align 8, !tbaa !1399, !dbg !2500
	store double  %359, double* %..inline.addr.740, align 8, !tbaa !1403, !dbg !2500
	br label %L..inline.13394, !dbg !2500
L..inline.13392:
	%360 = load i32, i32* %..inline.addr.550, align 4, !tbaa !1401, !dbg !2500

	%361 = sub i32  %360, 1, !dbg !2500
	store i32  %361, i32* %..inline.addr.550, align 4, !tbaa !1401, !dbg !2500
	%362 = sext i32  %360 to i64, !dbg !2500
	%363 = bitcast [20 x double]* %..inline.addr.620 to double*, !dbg !2500
	%364 = getelementptr double, double*  %363, i64  %362, !dbg !2500
	%365 = load double, double*  %364, align 8, !tbaa !1399, !dbg !2500
	store double  %365, double* %..inline.addr.740, align 8, !tbaa !1403, !dbg !2500
	br label %L..inline.13394
L..inline.13394:
	%366 = load double, double* %..inline.addr.740, align 8, !tbaa !1403, !dbg !2500
	store double  %366, double* %..inline.addr.600, align 8, !tbaa !1403, !dbg !2500
	%367 = load double, double* %..inline.addr.590, align 8, !tbaa !1403, !dbg !2500
	%368 = fadd double  %366,  %367, !dbg !2500
	store double  %368, double* %..inline.addr.590, align 8, !tbaa !1403, !dbg !2500
	%369 = load i32, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !2500

	%370 = add i32  %369, 1, !dbg !2500
	store i32  %370, i32* %..inline.addr.640, align 4, !tbaa !1401, !dbg !2500
	%371 = load i32, i32* %.x31746.addr, align 4, !tbaa !1401, !dbg !2500
	%372 = sub i32  %371, 1, !dbg !2500
	store i32  %372, i32* %.x31746.addr, align 4, !tbaa !1401, !dbg !2500
	%373 = icmp sgt i32  %372, 0, !dbg !2500
	br i1  %373, label %L..inline.13350, label %L..inline.13351, !dbg !2500
L..inline.13351:
	br label %L..inline.13344
L..inline.13344:
	%374 = load double, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2500
	%375 = call double @llvm.fabs.f64 (double  %374), !dbg !2500
	store double  %375, double* %..inline.addr.7, align 8, !tbaa !1403, !dbg !2500
	%376 = load double, double* %..inline.addr.600, align 8, !tbaa !1403, !dbg !2501
	%377 = call double @llvm.fabs.f64 (double  %376), !dbg !2501
	store double  %377, double* %..inline.addr.600, align 8, !tbaa !1403, !dbg !2501
	%378 = call double @llvm.maxnum.f64 (double  %375, double  %377), !dbg !2502
	%379 = fcmp uno double  %377,  %377, !dbg !2502
	%380 = select i1  %379, double  %375, double  %378, !dbg !2502
	%381 = load double, double* %..inline.addr.400, align 8, !tbaa !1403, !dbg !2503
	%382 = load double, double* %..inline.addr.590, align 8, !tbaa !1403, !dbg !2503
	%383 = fsub double  %381,  %382, !dbg !2503
	%384 = call double @llvm.fabs.f64 (double  %383), !dbg !2503
	%385 = call double @llvm.maxnum.f64 (double  %380, double  %384), !dbg !2504
	%386 = fcmp uno double  %384,  %384, !dbg !2504
	%387 = select i1  %386, double  %380, double  %385, !dbg !2504
	%388 = load double, double* %..inline.addr.4, align 8, !tbaa !1403, !dbg !2505
	%389 = fcmp uge double  %387,  %388, !dbg !2505
	br i1  %389, label %L..inline.13402, label %L.B0421, !dbg !2505
L.B0421:
	%390 = fadd double  %381,  %382, !dbg !2505
	%391 = fmul double  %390,  5.00000000000000000E-1, !dbg !2505
	store double  %391, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2505
	br label %L..inline.13403, !dbg !2506
L..inline.13402:
	br label %L..inline.13280
L..inline.13280:
	%392 = load i8*, i8** %.G0103.addr, align 8, !tbaa !1398, !dbg !2507
	%393 = getelementptr i8, i8*  %392, i64 18446744073709551608, !dbg !2507
	%394 = bitcast i8*  %393 to double*, !dbg !2507
	%395 = load double, double*  %394, align 8, !tbaa !1399, !dbg !2507
	%396 = bitcast i8*  %392 to double*, !dbg !2507
	store double  %395, double*  %396, align 8, !tbaa !1399, !dbg !2507
	%397 = load i8*, i8** %.G0101.addr, align 8, !tbaa !1398, !dbg !2507
	%398 = getelementptr i8, i8*  %397, i64 18446744073709551608, !dbg !2507
	%399 = bitcast i8*  %398 to double*, !dbg !2507
	%400 = load double, double*  %399, align 8, !tbaa !1399, !dbg !2507
	%401 = fmul double  %400,  2.50000000000000000E-1, !dbg !2507
	%402 = bitcast i8*  %397 to double*, !dbg !2507
	store double  %401, double*  %402, align 8, !tbaa !1399, !dbg !2507
	%403 = load i32, i32* %..inline.addr.8, align 4, !tbaa !1401, !dbg !2508

	%404 = add i32  %403, 1, !dbg !2509
	store i32  %404, i32* %..inline.addr.8, align 4, !tbaa !1401, !dbg !2509
	%405 = getelementptr i8, i8*  %392, i64 8, !dbg !2495
	store i8*  %405, i8** %.G0103.addr, align 8, !tbaa !1398, !dbg !2495
	%406 = getelementptr i8, i8*  %397, i64 8, !dbg !2495
	store i8*  %406, i8** %.G0101.addr, align 8, !tbaa !1398, !dbg !2495
	%407 = load i32, i32* %.x31745.addr, align 4, !tbaa !1401, !dbg !2495
	%408 = sub i32  %407, 1, !dbg !2495
	store i32  %408, i32* %.x31745.addr, align 4, !tbaa !1401, !dbg !2495
	%409 = icmp sgt i32  %408, 0, !dbg !2509
	br i1  %409, label %L..inline.13253, label %L.N0028, !dbg !2509
L.N0028:
	br label %L..inline.13403
L..inline.13403:
	%410 = load double, double* %..inline.addr.280, align 8, !tbaa !1403, !dbg !2510
	%411 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T2D_fix1 to i8*, !dbg !2490
	%412 = getelementptr i8, i8*  %411, i64 16, !dbg !2490
	%413 = bitcast %struct.T2D_fix1* %_T41506440_7661.addr to i8**, !dbg !2490
	store i8*  %412, i8**  %413, align 8, !tbaa !1398, !dbg !2490
	%414 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2511
	%415 = getelementptr i8, i8*  %414, i64 16, !dbg !2511
	store i8*  %415, i8**  %413, align 8, !tbaa !1398, !dbg !2511
	ret double  %410, !dbg !2490
L.B0422:
	br label %L.R0031, !dbg !2490
L_T41507328_7661:
	%416 = landingpad %astruct.dt64
	cleanup
	%417 = extractvalue %astruct.dt64  %416, 0, !dbg !2490
	store i8*  %417, i8** %__caught_object_address.addr, align 1, !tbaa !1398, !dbg !2490
	%418 = extractvalue %astruct.dt64  %416, 1, !dbg !2490
	store i32  %418, i32* %__catch_clause_number.addr, align 1, !tbaa !1399, !dbg !2490
	%419 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T2D_fix1 to i8*, !dbg !2490
	%420 = getelementptr i8, i8*  %419, i64 16, !dbg !2490
	%421 = bitcast %struct.T2D_fix1* %_T41506440_7661.addr to i8**, !dbg !2490
	store i8*  %420, i8**  %421, align 8, !tbaa !1398, !dbg !2490
	%422 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2511
	%423 = getelementptr i8, i8*  %422, i64 16, !dbg !2511
	store i8*  %423, i8**  %421, align 8, !tbaa !1398, !dbg !2511
	%424 = load i32, i32* %__catch_clause_number.addr, align 4, !tbaa !1399, !dbg !2490
	%425 = load i8*, i8** %__caught_object_address.addr, align 8, !tbaa !1398, !dbg !2490
	%426 = insertvalue %astruct.dt64 undef, i8*  %425, 0, !dbg !2490
	%427 = insertvalue %astruct.dt64  %426, i32  %424, 1, !dbg !2490
	resume %astruct.dt64  %427 , !dbg !2490
L.R0031:
	ret double 0.0
}
define linkonce_odr void @_ZN9Tinty_f2DD1Ev(%struct.Tinty_f2D* %_T41505256_7662.arg) #0 inlinehint !dbg !2515 {
L.entry:
	%_T41505256_7662.addr = alloca %struct.Tinty_f2D*, align 8

	call void @llvm.dbg.declare (metadata %struct.Tinty_f2D** %_T41505256_7662.addr, metadata !2527, metadata !1366), !dbg !2516
	store %struct.Tinty_f2D* %_T41505256_7662.arg, %struct.Tinty_f2D** %_T41505256_7662.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.Tinty_f2D** %_T41505256_7662.addr, metadata !2528, metadata !1366), !dbg !2516
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV9Tinty_f2D to i8*, !dbg !2529
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2529
	%2 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %_T41505256_7662.addr, align 8, !tbaa !1398, !dbg !2529
	%3 = bitcast %struct.Tinty_f2D*  %2 to i8**, !dbg !2529
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2529
	%4 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2529
	%5 = getelementptr i8, i8*  %4, i64 16, !dbg !2529
	store i8*  %5, i8**  %3, align 8, !tbaa !1398, !dbg !2529
	ret void, !dbg !2530
}
define linkonce_odr void @_ZN9Tinty_f2DD0Ev(%struct.Tinty_f2D* %_T41505256_7663.arg) #0 inlinehint !dbg !2532 {
L.entry:
	%_T41505256_7663.addr = alloca %struct.Tinty_f2D*, align 8
	%..inline.addr = alloca %struct.Tinty_f2D*, align 8

	store %struct.Tinty_f2D* %_T41505256_7663.arg, %struct.Tinty_f2D** %_T41505256_7663.addr, align 8, !tbaa !1398
	%0 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %_T41505256_7663.addr, align 8, !tbaa !1398, !dbg !2548
	%1 = bitcast %struct.Tinty_f2D*  %0 to i8*, !dbg !2548
	%2 = bitcast %struct.Tinty_f2D** %..inline.addr to i8**, !dbg !2548
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2548
	%3 = bitcast [5 x i32 (...)* (...)*]* @_ZTV9Tinty_f2D to i8*, !dbg !2548
	%4 = getelementptr i8, i8*  %3, i64 16, !dbg !2548
	%5 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %..inline.addr, align 8, !tbaa !1398, !dbg !2548
	%6 = bitcast %struct.Tinty_f2D*  %5 to i8**, !dbg !2548
	store i8*  %4, i8**  %6, align 8, !tbaa !1398, !dbg !2548
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2548
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2548
	store i8*  %8, i8**  %6, align 8, !tbaa !1398, !dbg !2548
	call void  @_ZdlPvm (i8*  %1, i64 40) nounwind, !dbg !2548
	ret void, !dbg !2548
}
define linkonce_odr void @_ZN9Tinty_f2DD2Ev(%struct.Tinty_f2D* %_T41505256_7664.arg) #0 inlinehint !dbg !2550 {
L.entry:
	%_T41505256_7664.addr = alloca %struct.Tinty_f2D*, align 8
	%..inline.addr = alloca %struct.Tinty_f2D*, align 8

	store %struct.Tinty_f2D* %_T41505256_7664.arg, %struct.Tinty_f2D** %_T41505256_7664.addr, align 8, !tbaa !1398
	%0 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %_T41505256_7664.addr, align 8, !tbaa !1398, !dbg !2566
	%1 = bitcast %struct.Tinty_f2D*  %0 to i8*, !dbg !2566
	%2 = bitcast %struct.Tinty_f2D** %..inline.addr to i8**, !dbg !2566
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2566
	%3 = bitcast [5 x i32 (...)* (...)*]* @_ZTV9Tinty_f2D to i8*, !dbg !2566
	%4 = getelementptr i8, i8*  %3, i64 16, !dbg !2566
	%5 = load %struct.Tinty_f2D*, %struct.Tinty_f2D** %..inline.addr, align 8, !tbaa !1398, !dbg !2566
	%6 = bitcast %struct.Tinty_f2D*  %5 to i8**, !dbg !2566
	store i8*  %4, i8**  %6, align 8, !tbaa !1398, !dbg !2566
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2566
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2566
	store i8*  %8, i8**  %6, align 8, !tbaa !1398, !dbg !2566
	ret void, !dbg !2566
}
define linkonce_odr void @_ZN10Tintxy_f3DC1ERK11T3DFunctionddddd(%struct.Tintxy_f3D* %_T41505256_7665.arg, %struct.T3DFunction* %f1.arg, double %xmin.arg, double %xmax.arg, double %ymin.arg, double %ymax.arg, double %absacc.arg) #0 inlinehint !dbg !2570 {
L.entry:
	%_T41505256_7665.addr = alloca %struct.Tintxy_f3D*, align 8
	%f1.addr = alloca %struct.T3DFunction*, align 8
	%xmin.addr = alloca double, align 8
	%xmax.addr = alloca double, align 8
	%ymin.addr = alloca double, align 8
	%ymax.addr = alloca double, align 8
	%absacc.addr = alloca double, align 8

	call void @llvm.dbg.declare (metadata %struct.Tintxy_f3D** %_T41505256_7665.addr, metadata !2582, metadata !1366), !dbg !2571
	store %struct.Tintxy_f3D* %_T41505256_7665.arg, %struct.Tintxy_f3D** %_T41505256_7665.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.Tintxy_f3D** %_T41505256_7665.addr, metadata !2583, metadata !1366), !dbg !2571
	call void @llvm.dbg.declare (metadata %struct.T3DFunction** %f1.addr, metadata !2584, metadata !1366), !dbg !2571
	store %struct.T3DFunction* %f1.arg, %struct.T3DFunction** %f1.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.T3DFunction** %f1.addr, metadata !2585, metadata !1366), !dbg !2571
	call void @llvm.dbg.declare (metadata double* %xmin.addr, metadata !2586, metadata !1366), !dbg !2571
	store double %xmin.arg, double* %xmin.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %xmin.addr, metadata !2587, metadata !1366), !dbg !2571
	call void @llvm.dbg.declare (metadata double* %xmax.addr, metadata !2588, metadata !1366), !dbg !2571
	store double %xmax.arg, double* %xmax.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %xmax.addr, metadata !2589, metadata !1366), !dbg !2571
	call void @llvm.dbg.declare (metadata double* %ymin.addr, metadata !2590, metadata !1366), !dbg !2571
	store double %ymin.arg, double* %ymin.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %ymin.addr, metadata !2591, metadata !1366), !dbg !2571
	call void @llvm.dbg.declare (metadata double* %ymax.addr, metadata !2592, metadata !1366), !dbg !2571
	store double %ymax.arg, double* %ymax.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %ymax.addr, metadata !2593, metadata !1366), !dbg !2571
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !2594, metadata !1366), !dbg !2571
	store double %absacc.arg, double* %absacc.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %absacc.addr, metadata !2595, metadata !1366), !dbg !2571
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2596
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2596
	%2 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %_T41505256_7665.addr, align 8, !tbaa !1398, !dbg !2596
	%3 = bitcast %struct.Tintxy_f3D*  %2 to i8**, !dbg !2596
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2596
	%4 = bitcast [5 x i32 (...)* (...)*]* @_ZTV10Tintxy_f3D to i8*, !dbg !2596
	%5 = getelementptr i8, i8*  %4, i64 16, !dbg !2596
	store i8*  %5, i8**  %3, align 8, !tbaa !1398, !dbg !2596
	%6 = load %struct.T3DFunction*, %struct.T3DFunction** %f1.addr, align 8, !tbaa !1398, !dbg !2596
	%7 = bitcast %struct.T3DFunction*  %6 to i8*, !dbg !2596
	%8 = bitcast %struct.Tintxy_f3D*  %2 to i8*, !dbg !2596
	%9 = getelementptr i8, i8*  %8, i64 8, !dbg !2596
	%10 = bitcast i8*  %9 to i8**, !dbg !2596
	store i8*  %7, i8**  %10, align 8, !tbaa !1398, !dbg !2596
	%11 = load double, double* %xmin.addr, align 8, !tbaa !1403, !dbg !2596
	%12 = getelementptr i8, i8*  %8, i64 16, !dbg !2596
	%13 = bitcast i8*  %12 to double*, !dbg !2596
	store double  %11, double*  %13, align 8, !tbaa !1399, !dbg !2596
	%14 = load double, double* %xmax.addr, align 8, !tbaa !1403, !dbg !2596
	%15 = getelementptr i8, i8*  %8, i64 24, !dbg !2596
	%16 = bitcast i8*  %15 to double*, !dbg !2596
	store double  %14, double*  %16, align 8, !tbaa !1399, !dbg !2596
	%17 = load double, double* %ymin.addr, align 8, !tbaa !1403, !dbg !2596
	%18 = getelementptr i8, i8*  %8, i64 32, !dbg !2596
	%19 = bitcast i8*  %18 to double*, !dbg !2596
	store double  %17, double*  %19, align 8, !tbaa !1399, !dbg !2596
	%20 = load double, double* %ymax.addr, align 8, !tbaa !1403, !dbg !2596
	%21 = getelementptr i8, i8*  %8, i64 40, !dbg !2596
	%22 = bitcast i8*  %21 to double*, !dbg !2596
	store double  %20, double*  %22, align 8, !tbaa !1399, !dbg !2596
	%23 = load double, double* %absacc.addr, align 8, !tbaa !1403, !dbg !2596
	%24 = getelementptr i8, i8*  %8, i64 48, !dbg !2596
	%25 = bitcast i8*  %24 to double*, !dbg !2596
	store double  %23, double*  %25, align 8, !tbaa !1399, !dbg !2596
	ret void, !dbg !2596
}
define linkonce_odr void @_ZN10Tintxy_f3DC2ERK11T3DFunctionddddd(%struct.Tintxy_f3D* %_T41505256_7666.arg, %struct.T3DFunction* %_T41505552_7666.arg, double %_T41505848_7666.arg, double %_T41506144_7666.arg, double %_T41506440_7666.arg, double %_T41506736_7666.arg, double %_T41507032_7666.arg) #0 inlinehint !dbg !2598 {
L.entry:
	%_T41505256_7666.addr = alloca %struct.Tintxy_f3D*, align 8
	%_T41505552_7666.addr = alloca %struct.T3DFunction*, align 8
	%_T41505848_7666.addr = alloca double, align 8
	%_T41506144_7666.addr = alloca double, align 8
	%_T41506440_7666.addr = alloca double, align 8
	%_T41506736_7666.addr = alloca double, align 8
	%_T41507032_7666.addr = alloca double, align 8
	%..inline.addr = alloca %struct.Tintxy_f3D*, align 8
	%..inline.addr.1 = alloca %struct.T3DFunction*, align 8
	%..inline.addr.2 = alloca double, align 8
	%..inline.addr.3 = alloca double, align 8
	%..inline.addr.4 = alloca double, align 8
	%..inline.addr.5 = alloca double, align 8
	%..inline.addr.6 = alloca double, align 8

	store %struct.Tintxy_f3D* %_T41505256_7666.arg, %struct.Tintxy_f3D** %_T41505256_7666.addr, align 8, !tbaa !1398
	store %struct.T3DFunction* %_T41505552_7666.arg, %struct.T3DFunction** %_T41505552_7666.addr, align 8, !tbaa !1398
	store double %_T41505848_7666.arg, double* %_T41505848_7666.addr, align 8, !tbaa !1399
	store double %_T41506144_7666.arg, double* %_T41506144_7666.addr, align 8, !tbaa !1399
	store double %_T41506440_7666.arg, double* %_T41506440_7666.addr, align 8, !tbaa !1399
	store double %_T41506736_7666.arg, double* %_T41506736_7666.addr, align 8, !tbaa !1399
	store double %_T41507032_7666.arg, double* %_T41507032_7666.addr, align 8, !tbaa !1399
	%0 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %_T41505256_7666.addr, align 8, !tbaa !1398, !dbg !2614
	%1 = bitcast %struct.Tintxy_f3D*  %0 to i8*, !dbg !2614
	%2 = bitcast %struct.Tintxy_f3D** %..inline.addr to i8**, !dbg !2614
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2614
	%3 = load %struct.T3DFunction*, %struct.T3DFunction** %_T41505552_7666.addr, align 8, !tbaa !1398, !dbg !2614
	%4 = bitcast %struct.T3DFunction*  %3 to i8*, !dbg !2614
	%5 = bitcast %struct.T3DFunction** %..inline.addr.1 to i8**, !dbg !2614
	store i8*  %4, i8**  %5, align 8, !tbaa !1398, !dbg !2614
	%6 = load double, double* %_T41505848_7666.addr, align 8, !tbaa !1403, !dbg !2614
	%7 = load double, double* %_T41506144_7666.addr, align 8, !tbaa !1403, !dbg !2614
	%8 = load double, double* %_T41506440_7666.addr, align 8, !tbaa !1403, !dbg !2614
	%9 = load double, double* %_T41506736_7666.addr, align 8, !tbaa !1403, !dbg !2614
	%10 = load double, double* %_T41507032_7666.addr, align 8, !tbaa !1403, !dbg !2614
	%11 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2614
	%12 = getelementptr i8, i8*  %11, i64 16, !dbg !2614
	%13 = bitcast %struct.Tintxy_f3D*  %0 to i8**, !dbg !2614
	store i8*  %12, i8**  %13, align 8, !tbaa !1398, !dbg !2614
	%14 = bitcast [5 x i32 (...)* (...)*]* @_ZTV10Tintxy_f3D to i8*, !dbg !2614
	%15 = getelementptr i8, i8*  %14, i64 16, !dbg !2614
	%16 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %..inline.addr, align 8, !tbaa !1398, !dbg !2614
	%17 = bitcast %struct.Tintxy_f3D*  %16 to i8**, !dbg !2614
	store i8*  %15, i8**  %17, align 8, !tbaa !1398, !dbg !2614
	%18 = load %struct.T3DFunction*, %struct.T3DFunction** %..inline.addr.1, align 8, !tbaa !1398, !dbg !2614
	%19 = bitcast %struct.T3DFunction*  %18 to i8*, !dbg !2614
	%20 = bitcast %struct.Tintxy_f3D*  %16 to i8*, !dbg !2614
	%21 = getelementptr i8, i8*  %20, i64 8, !dbg !2614
	%22 = bitcast i8*  %21 to i8**, !dbg !2614
	store i8*  %19, i8**  %22, align 8, !tbaa !1398, !dbg !2614
	%23 = getelementptr i8, i8*  %20, i64 16, !dbg !2614
	%24 = bitcast i8*  %23 to double*, !dbg !2614
	store double  %6, double*  %24, align 8, !tbaa !1399, !dbg !2614
	%25 = getelementptr i8, i8*  %20, i64 24, !dbg !2614
	%26 = bitcast i8*  %25 to double*, !dbg !2614
	store double  %7, double*  %26, align 8, !tbaa !1399, !dbg !2614
	%27 = getelementptr i8, i8*  %20, i64 32, !dbg !2614
	%28 = bitcast i8*  %27 to double*, !dbg !2614
	store double  %8, double*  %28, align 8, !tbaa !1399, !dbg !2614
	%29 = getelementptr i8, i8*  %20, i64 40, !dbg !2614
	%30 = bitcast i8*  %29 to double*, !dbg !2614
	store double  %9, double*  %30, align 8, !tbaa !1399, !dbg !2614
	%31 = getelementptr i8, i8*  %20, i64 48, !dbg !2614
	%32 = bitcast i8*  %31 to double*, !dbg !2614
	store double  %10, double*  %32, align 8, !tbaa !1399, !dbg !2614
	ret void, !dbg !2614
}
define linkonce_odr double @_ZNK10Tintxy_f3D4callEd(%struct.Tintxy_f3D* %_T41505256_7667.arg, double %z.arg) #0 inlinehint personality i8* bitcast (i32 (...)* @__gxx_personality_v0 to i8*) !dbg !2616 {
L.entry:
	%_T41505256_7667.addr = alloca %struct.Tintxy_f3D*, align 8
	%z.addr = alloca double, align 8
	%..inline.addr = alloca %struct.T3DFunction*, align 8
	%..inline.addr.1 = alloca double, align 8
	%_T41506440_7667.addr = alloca %struct.T3D_fix3, align 8
	%.Q0025.addr = alloca double, align 8
	%_T41506736_7667.addr = alloca double, align 8
	%__caught_object_address.addr = alloca i8*, align 8
	%__catch_clause_number.addr = alloca i32, align 4

	call void @llvm.dbg.declare (metadata %struct.Tintxy_f3D** %_T41505256_7667.addr, metadata !2656, metadata !1366), !dbg !2617
	store %struct.Tintxy_f3D* %_T41505256_7667.arg, %struct.Tintxy_f3D** %_T41505256_7667.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.Tintxy_f3D** %_T41505256_7667.addr, metadata !2657, metadata !1366), !dbg !2617
	call void @llvm.dbg.declare (metadata double* %z.addr, metadata !2658, metadata !1366), !dbg !2617
	store double %z.arg, double* %z.addr, align 8, !tbaa !1399
	call void @llvm.dbg.declare (metadata double* %z.addr, metadata !2659, metadata !1366), !dbg !2617
	%0 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %_T41505256_7667.addr, align 8, !tbaa !1398, !dbg !2660
	%1 = bitcast %struct.Tintxy_f3D*  %0 to i8*, !dbg !2660
	%2 = getelementptr i8, i8*  %1, i64 8, !dbg !2660
	%3 = bitcast i8*  %2 to i8**, !dbg !2660
	%4 = load i8*, i8**  %3, align 8, !tbaa !1398, !dbg !2660
	%5 = bitcast %struct.T3DFunction** %..inline.addr to i8**, !dbg !2660
	store i8*  %4, i8**  %5, align 8, !tbaa !1398, !dbg !2660
	%6 = load double, double* %z.addr, align 8, !tbaa !1403, !dbg !2660
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2661
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2661
	%9 = bitcast %struct.T3D_fix3* %_T41506440_7667.addr to i8**, !dbg !2661
	store i8*  %8, i8**  %9, align 8, !tbaa !1398, !dbg !2661
	%10 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T3D_fix3 to i8*, !dbg !2661
	%11 = getelementptr i8, i8*  %10, i64 16, !dbg !2661
	store i8*  %11, i8**  %9, align 8, !tbaa !1398, !dbg !2661
	%12 = load %struct.T3DFunction*, %struct.T3DFunction** %..inline.addr, align 8, !tbaa !1398, !dbg !2661
	%13 = bitcast %struct.T3DFunction*  %12 to i8*, !dbg !2661
	%14 = bitcast %struct.T3D_fix3* %_T41506440_7667.addr to i8*, !dbg !2661
	%15 = getelementptr i8, i8*  %14, i64 8, !dbg !2661
	%16 = bitcast i8*  %15 to i8**, !dbg !2661
	store i8*  %13, i8**  %16, align 8, !tbaa !1398, !dbg !2661
	%17 = getelementptr i8, i8*  %14, i64 16, !dbg !2661
	%18 = bitcast i8*  %17 to double*, !dbg !2661
	store double  %6, double*  %18, align 8, !tbaa !1399, !dbg !2661
	%19 = bitcast %struct.T3D_fix3* %_T41506440_7667.addr to %struct.T2DFunction*, !dbg !2660
	%20 = getelementptr i8, i8*  %1, i64 16, !dbg !2660
	%21 = bitcast i8*  %20 to double*, !dbg !2660
	%22 = load double, double*  %21, align 8, !tbaa !1399, !dbg !2660
	%23 = getelementptr i8, i8*  %1, i64 24, !dbg !2660
	%24 = bitcast i8*  %23 to double*, !dbg !2660
	%25 = load double, double*  %24, align 8, !tbaa !1399, !dbg !2660
	%26 = getelementptr i8, i8*  %1, i64 32, !dbg !2660
	%27 = bitcast i8*  %26 to double*, !dbg !2660
	%28 = load double, double*  %27, align 8, !tbaa !1399, !dbg !2660
	%29 = getelementptr i8, i8*  %1, i64 40, !dbg !2660
	%30 = bitcast i8*  %29 to double*, !dbg !2660
	%31 = load double, double*  %30, align 8, !tbaa !1399, !dbg !2660
	%32 = getelementptr i8, i8*  %1, i64 48, !dbg !2660
	%33 = bitcast i8*  %32 to double*, !dbg !2660
	%34 = load double, double*  %33, align 8, !tbaa !1399, !dbg !2660
	%35 = invoke double  @_Z7RombergRK11T2DFunctionddddd (%struct.T2DFunction*  %19, double  %22, double  %25, double  %28, double  %31, double  %34)
		to label %L.B0439
		unwind label %L_T41507328_7667, !dbg !2660
L.B0439:
	store double  %35, double* %.Q0025.addr, align 8, !tbaa !1403, !dbg !2660
	%36 = getelementptr i8, i8*  %10, i64 16, !dbg !2660
	store i8*  %36, i8**  %9, align 8, !tbaa !1398, !dbg !2660
	%37 = getelementptr i8, i8*  %7, i64 16, !dbg !2662
	store i8*  %37, i8**  %9, align 8, !tbaa !1398, !dbg !2662
	ret double  %35, !dbg !2660
L.B0440:
	br label %L.R0037, !dbg !2660
L_T41507328_7667:
	%38 = landingpad %astruct.dt64
	cleanup
	%39 = extractvalue %astruct.dt64  %38, 0, !dbg !2660
	store i8*  %39, i8** %__caught_object_address.addr, align 1, !tbaa !1398, !dbg !2660
	%40 = extractvalue %astruct.dt64  %38, 1, !dbg !2660
	store i32  %40, i32* %__catch_clause_number.addr, align 1, !tbaa !1399, !dbg !2660
	%41 = bitcast [5 x i32 (...)* (...)*]* @_ZTV8T3D_fix3 to i8*, !dbg !2660
	%42 = getelementptr i8, i8*  %41, i64 16, !dbg !2660
	%43 = bitcast %struct.T3D_fix3* %_T41506440_7667.addr to i8**, !dbg !2660
	store i8*  %42, i8**  %43, align 8, !tbaa !1398, !dbg !2660
	%44 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T2DFunction to i8*, !dbg !2662
	%45 = getelementptr i8, i8*  %44, i64 16, !dbg !2662
	store i8*  %45, i8**  %43, align 8, !tbaa !1398, !dbg !2662
	%46 = load i32, i32* %__catch_clause_number.addr, align 4, !tbaa !1399, !dbg !2660
	%47 = load i8*, i8** %__caught_object_address.addr, align 8, !tbaa !1398, !dbg !2660
	%48 = insertvalue %astruct.dt64 undef, i8*  %47, 0, !dbg !2660
	%49 = insertvalue %astruct.dt64  %48, i32  %46, 1, !dbg !2660
	resume %astruct.dt64  %49 , !dbg !2660
L.R0037:
	ret double 0.0
}
define linkonce_odr void @_ZN10Tintxy_f3DD1Ev(%struct.Tintxy_f3D* %_T41505256_7668.arg) #0 inlinehint !dbg !2666 {
L.entry:
	%_T41505256_7668.addr = alloca %struct.Tintxy_f3D*, align 8

	call void @llvm.dbg.declare (metadata %struct.Tintxy_f3D** %_T41505256_7668.addr, metadata !2678, metadata !1366), !dbg !2667
	store %struct.Tintxy_f3D* %_T41505256_7668.arg, %struct.Tintxy_f3D** %_T41505256_7668.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct.Tintxy_f3D** %_T41505256_7668.addr, metadata !2679, metadata !1366), !dbg !2667
	%0 = bitcast [5 x i32 (...)* (...)*]* @_ZTV10Tintxy_f3D to i8*, !dbg !2680
	%1 = getelementptr i8, i8*  %0, i64 16, !dbg !2680
	%2 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %_T41505256_7668.addr, align 8, !tbaa !1398, !dbg !2680
	%3 = bitcast %struct.Tintxy_f3D*  %2 to i8**, !dbg !2680
	store i8*  %1, i8**  %3, align 8, !tbaa !1398, !dbg !2680
	%4 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2680
	%5 = getelementptr i8, i8*  %4, i64 16, !dbg !2680
	store i8*  %5, i8**  %3, align 8, !tbaa !1398, !dbg !2680
	ret void, !dbg !2681
}
define linkonce_odr void @_ZN10Tintxy_f3DD0Ev(%struct.Tintxy_f3D* %_T41505256_7669.arg) #0 inlinehint !dbg !2683 {
L.entry:
	%_T41505256_7669.addr = alloca %struct.Tintxy_f3D*, align 8
	%..inline.addr = alloca %struct.Tintxy_f3D*, align 8

	store %struct.Tintxy_f3D* %_T41505256_7669.arg, %struct.Tintxy_f3D** %_T41505256_7669.addr, align 8, !tbaa !1398
	%0 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %_T41505256_7669.addr, align 8, !tbaa !1398, !dbg !2699
	%1 = bitcast %struct.Tintxy_f3D*  %0 to i8*, !dbg !2699
	%2 = bitcast %struct.Tintxy_f3D** %..inline.addr to i8**, !dbg !2699
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2699
	%3 = bitcast [5 x i32 (...)* (...)*]* @_ZTV10Tintxy_f3D to i8*, !dbg !2699
	%4 = getelementptr i8, i8*  %3, i64 16, !dbg !2699
	%5 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %..inline.addr, align 8, !tbaa !1398, !dbg !2699
	%6 = bitcast %struct.Tintxy_f3D*  %5 to i8**, !dbg !2699
	store i8*  %4, i8**  %6, align 8, !tbaa !1398, !dbg !2699
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2699
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2699
	store i8*  %8, i8**  %6, align 8, !tbaa !1398, !dbg !2699
	call void  @_ZdlPvm (i8*  %1, i64 56) nounwind, !dbg !2699
	ret void, !dbg !2699
}
define linkonce_odr void @_ZN10Tintxy_f3DD2Ev(%struct.Tintxy_f3D* %_T41505256_7670.arg) #0 inlinehint !dbg !2701 {
L.entry:
	%_T41505256_7670.addr = alloca %struct.Tintxy_f3D*, align 8
	%..inline.addr = alloca %struct.Tintxy_f3D*, align 8

	store %struct.Tintxy_f3D* %_T41505256_7670.arg, %struct.Tintxy_f3D** %_T41505256_7670.addr, align 8, !tbaa !1398
	%0 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %_T41505256_7670.addr, align 8, !tbaa !1398, !dbg !2717
	%1 = bitcast %struct.Tintxy_f3D*  %0 to i8*, !dbg !2717
	%2 = bitcast %struct.Tintxy_f3D** %..inline.addr to i8**, !dbg !2717
	store i8*  %1, i8**  %2, align 8, !tbaa !1398, !dbg !2717
	%3 = bitcast [5 x i32 (...)* (...)*]* @_ZTV10Tintxy_f3D to i8*, !dbg !2717
	%4 = getelementptr i8, i8*  %3, i64 16, !dbg !2717
	%5 = load %struct.Tintxy_f3D*, %struct.Tintxy_f3D** %..inline.addr, align 8, !tbaa !1398, !dbg !2717
	%6 = bitcast %struct.Tintxy_f3D*  %5 to i8**, !dbg !2717
	store i8*  %4, i8**  %6, align 8, !tbaa !1398, !dbg !2717
	%7 = bitcast [5 x i32 (...)* (...)*]* @_ZTV11T1DFunction to i8*, !dbg !2717
	%8 = getelementptr i8, i8*  %7, i64 16, !dbg !2717
	store i8*  %8, i8**  %6, align 8, !tbaa !1398, !dbg !2717
	ret void, !dbg !2717
}
define linkonce_odr i64 @_ZNSt11char_traitsIcE6lengthEPKc(i8* %__s.arg) #0 inlinehint !dbg !2723 {
L.entry:
	%__s.addr = alloca i8*, align 8

	call void @llvm.dbg.declare (metadata i8** %__s.addr, metadata !2727, metadata !1366), !dbg !2724
	store i8* %__s.arg, i8** %__s.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata i8** %__s.addr, metadata !2728, metadata !1366), !dbg !2724
	%0 = load i8*, i8** %__s.addr, align 8, !tbaa !1398, !dbg !2729
	%1 = call i64  @strlen (i8*  %0) nounwind, !dbg !2729
	ret i64  %1, !dbg !2729
}
define linkonce_odr signext i32 @_ZNKSt9basic_iosIcSt11char_traitsIcEE7rdstateEv(%struct._ZSt9basic_iosIcSt11char_traitsIcEE* %_T41505256_7676.arg) #0 inlinehint !dbg !2737 {
L.entry:
	%_T41505256_7676.addr = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8

	call void @llvm.dbg.declare (metadata %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %_T41505256_7676.addr, metadata !2741, metadata !1366), !dbg !2738
	store %struct._ZSt9basic_iosIcSt11char_traitsIcEE* %_T41505256_7676.arg, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %_T41505256_7676.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %_T41505256_7676.addr, metadata !2742, metadata !1366), !dbg !2738
	%0 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %_T41505256_7676.addr, align 8, !tbaa !1398, !dbg !2743
	%1 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %0 to i8*, !dbg !2743
	%2 = getelementptr i8, i8*  %1, i64 32, !dbg !2743
	%3 = bitcast i8*  %2 to i32*, !dbg !2743
	%4 = load i32, i32*  %3, align 4, !tbaa !1399, !dbg !2743
	ret i32  %4, !dbg !2743
}
define linkonce_odr void @_ZNSt9basic_iosIcSt11char_traitsIcEE8setstateESt12_Ios_Iostate(%struct._ZSt9basic_iosIcSt11char_traitsIcEE* %_T41505256_7680.arg, i32 signext %__state.arg) #0 inlinehint !dbg !2747 {
L.entry:
	%_T41505256_7680.addr = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%__state.addr = alloca i32, align 4
	%..inline.addr = alloca i32, align 4

	call void @llvm.dbg.declare (metadata %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %_T41505256_7680.addr, metadata !2759, metadata !1366), !dbg !2748
	store %struct._ZSt9basic_iosIcSt11char_traitsIcEE* %_T41505256_7680.arg, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %_T41505256_7680.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %_T41505256_7680.addr, metadata !2760, metadata !1366), !dbg !2748
	call void @llvm.dbg.declare (metadata i32* %__state.addr, metadata !2761, metadata !1366), !dbg !2748
	store i32 %__state.arg, i32* %__state.addr, align 4, !tbaa !1399
	call void @llvm.dbg.declare (metadata i32* %__state.addr, metadata !2762, metadata !1366), !dbg !2748
	%0 = load i32, i32* %__state.addr, align 4, !tbaa !1399, !dbg !2763
	%1 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %_T41505256_7680.addr, align 8, !tbaa !1398, !dbg !2763
	%2 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %1 to i8*, !dbg !2763
	%3 = getelementptr i8, i8*  %2, i64 32, !dbg !2763
	%4 = bitcast i8*  %3 to i32*, !dbg !2763
	%5 = load i32, i32*  %4, align 4, !tbaa !1399, !dbg !2763
	%6 = or i32  %0,  %5, !dbg !2763
	call void  @_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate (%struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %1, i32  %6), !dbg !2764
	ret void, !dbg !2763
}
define linkonce_odr signext i32 @_ZStorSt12_Ios_IostateS_(i32 signext %__a.arg, i32 signext %__b.arg) #0 inlinehint !dbg !2770 {
L.entry:
	%__a.addr = alloca i32, align 4
	%__b.addr = alloca i32, align 4

	call void @llvm.dbg.declare (metadata i32* %__a.addr, metadata !2774, metadata !1366), !dbg !2771
	store i32 %__a.arg, i32* %__a.addr, align 4, !tbaa !1399
	call void @llvm.dbg.declare (metadata i32* %__a.addr, metadata !2775, metadata !1366), !dbg !2771
	call void @llvm.dbg.declare (metadata i32* %__b.addr, metadata !2776, metadata !1366), !dbg !2771
	store i32 %__b.arg, i32* %__b.addr, align 4, !tbaa !1399
	call void @llvm.dbg.declare (metadata i32* %__b.addr, metadata !2777, metadata !1366), !dbg !2771
	%0 = load i32, i32* %__b.addr, align 4, !tbaa !1399, !dbg !2778
	%1 = load i32, i32* %__a.addr, align 4, !tbaa !1399, !dbg !2778
	%2 = or i32  %0,  %1, !dbg !2778
	ret i32  %2, !dbg !2778
}
define linkonce_odr %struct._ZSo* @_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc(%struct._ZSo* %__out.arg, i8* %__s.arg) #0 inlinehint !dbg !2789 {
L.entry:
	%__out.addr = alloca %struct._ZSo*, align 8
	%__s.addr = alloca i8*, align 8
	%..inline.addr = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.1 = alloca %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, align 8
	%..inline.addr.2 = alloca i32, align 4
	%..inline.addr.3 = alloca i8*, align 8
	%..inline.addr.4 = alloca i64, align 8
	%.Q0026.addr = alloca %struct._ZSo*, align 8

	call void @llvm.dbg.declare (metadata %struct._ZSo** %__out.addr, metadata !2809, metadata !1366), !dbg !2790
	store %struct._ZSo* %__out.arg, %struct._ZSo** %__out.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata %struct._ZSo** %__out.addr, metadata !2810, metadata !1366), !dbg !2790
	call void @llvm.dbg.declare (metadata i8** %__s.addr, metadata !2811, metadata !1366), !dbg !2790
	store i8* %__s.arg, i8** %__s.addr, align 8, !tbaa !1398
	call void @llvm.dbg.declare (metadata i8** %__s.addr, metadata !2812, metadata !1366), !dbg !2790
	%0 = load i8*, i8** %__s.addr, align 8, !tbaa !1398, !dbg !2813
	%1 = icmp ne i8*  %0,  null, !dbg !2813
	br i1  %1, label %L.B0048, label %L.B0455, !dbg !2813
L.B0455:
	%2 = load %struct._ZSo*, %struct._ZSo** %__out.addr, align 8, !tbaa !1398, !dbg !2814
	%3 = bitcast %struct._ZSo*  %2 to i8*, !dbg !2814
	%4 = bitcast %struct._ZSo*  %2 to i8**, !dbg !2814
	%5 = load i8*, i8**  %4, align 8, !tbaa !1398, !dbg !2814
	%6 = getelementptr i8, i8*  %5, i64 18446744073709551592, !dbg !2814
	%7 = bitcast i8*  %6 to i64*, !dbg !2814
	%8 = load i64, i64*  %7, align 8, !tbaa !1399, !dbg !2814
	%9 = getelementptr i8, i8*  %3, i64  %8, !dbg !2814
	%10 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr to i8**, !dbg !2814
	store i8*  %9, i8**  %10, align 8, !tbaa !1398, !dbg !2814
	%11 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr, align 8, !tbaa !1398, !dbg !2815
	%12 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %11 to i8*, !dbg !2815
	%13 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.1 to i8**, !dbg !2815
	store i8*  %12, i8**  %13, align 8, !tbaa !1398, !dbg !2815
	%14 = load %struct._ZSt9basic_iosIcSt11char_traitsIcEE*, %struct._ZSt9basic_iosIcSt11char_traitsIcEE** %..inline.addr.1, align 8, !tbaa !1398, !dbg !2815
	%15 = bitcast %struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %14 to i8*, !dbg !2815
	%16 = getelementptr i8, i8*  %15, i64 32, !dbg !2815
	%17 = bitcast i8*  %16 to i32*, !dbg !2815
	%18 = load i32, i32*  %17, align 4, !tbaa !1399, !dbg !2815
	%19 = or i32  %18, 1, !dbg !2815
	call void  @_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate (%struct._ZSt9basic_iosIcSt11char_traitsIcEE*  %11, i32  %19), !dbg !2815
	br label %L.B0049, !dbg !2814
L.B0048:
	%20 = load i8*, i8** %__s.addr, align 8, !tbaa !1398, !dbg !2816
	%21 = call i64  @strlen (i8*  %20) nounwind, !dbg !2817
	%22 = load %struct._ZSo*, %struct._ZSo** %__out.addr, align 8, !tbaa !1398, !dbg !2818
	%23 = load i8*, i8** %__s.addr, align 8, !tbaa !1398, !dbg !2818
	%24 = call %struct._ZSo*  @_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l (%struct._ZSo*  %22, i8*  %23, i64  %21), !dbg !2818
	%25 = bitcast %struct._ZSo*  %24 to i8*, !dbg !2818
	%26 = bitcast %struct._ZSo** %.Q0026.addr to i8**, !dbg !2818
	store i8*  %25, i8**  %26, align 8, !tbaa !1398, !dbg !2818
	br label %L.B0049
L.B0049:
	%27 = load %struct._ZSo*, %struct._ZSo** %__out.addr, align 8, !tbaa !1398, !dbg !2819
	ret %struct._ZSo*  %27, !dbg !2819
}

%struct._ZNSt8ios_base4InitE = type <{ [1 x i8]}> 

define void @__sti___25_backgroundfield_quadr_cpp_6e0a9365() #0 inlinehint !dbg !2822 {
L.entry:

	%0 = load i32, i32* @__I___25_backgroundfield_quadr_cpp_6e0a9365, align 4, !tbaa !1401, !dbg !2826
	%1 = icmp eq i32  %0, 1, !dbg !2826
	br i1  %1, label %L.B0050, label %L.B0456, !dbg !2826
L.B0456:
	store i32 1, i32* @__I___25_backgroundfield_quadr_cpp_6e0a9365, align 4, !tbaa !1401, !dbg !2826
	call void  @_ZNSt8ios_base4InitC1Ev (%struct._ZNSt8ios_base4InitE* @_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE) nounwind, !dbg !2826
	%2 = bitcast void (%struct._ZNSt8ios_base4InitE*)* @_ZNSt8ios_base4InitD1Ev to void (i8*)*, !dbg !2826
	%3 = bitcast %struct._ZNSt8ios_base4InitE* @_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE to i8*, !dbg !2826
	%4 = bitcast i8** @__dso_handle to i8*, !dbg !2826
	%5 = call i32  @__cxa_atexit (void (i8*)*  %2, i8*  %3, i8*  %4) nounwind, !dbg !2826
	br label %L.B0050
L.B0050:
	ret void, !dbg !2826
}
@_ZTV11T1DFunction = weak unnamed_addr global [5 x i32 (...)* (...)*]  [ i32 (...)* (...)*  null, i32 (...)* (...)*  bitcast (%struct.__class_type_info* @_ZTI11T1DFunction to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void ()* @__cxa_pure_virtual to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.T1DFunction*)* @_ZN11T1DFunctionD1Ev to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.T1DFunction*)* @_ZN11T1DFunctionD0Ev to i32 (...)* (...)*) ], align 16, !dbg !1871

%struct.__EDG_type_info = type <{ i32 (...)* (...)*, i8*}> 
%struct.__class_type_info = type <{ %struct.__EDG_type_info}> 

@_ZTV11T2DFunction = weak unnamed_addr global [5 x i32 (...)* (...)*]  [ i32 (...)* (...)*  null, i32 (...)* (...)*  bitcast (%struct.__class_type_info* @_ZTI11T2DFunction to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void ()* @__cxa_pure_virtual to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.T2DFunction*)* @_ZN11T2DFunctionD1Ev to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.T2DFunction*)* @_ZN11T2DFunctionD0Ev to i32 (...)* (...)*) ], align 16, !dbg !2096
@_ZTV8T2D_fix1 = weak unnamed_addr global [5 x i32 (...)* (...)*]  [ i32 (...)* (...)*  null, i32 (...)* (...)*  bitcast (%struct.__si_class_type_info* @_ZTI8T2D_fix1 to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (double (%struct.T2D_fix1*, double)* @_ZNK8T2D_fix14callEd to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.T2D_fix1*)* @_ZN8T2D_fix1D1Ev to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.T2D_fix1*)* @_ZN8T2D_fix1D0Ev to i32 (...)* (...)*) ], align 16, !dbg !2158

%struct.__si_class_type_info = type <{ %struct.__class_type_info, %struct.__class_type_info*}> 

@_ZTV8T3D_fix3 = weak unnamed_addr global [5 x i32 (...)* (...)*]  [ i32 (...)* (...)*  null, i32 (...)* (...)*  bitcast (%struct.__si_class_type_info* @_ZTI8T3D_fix3 to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (double (%struct.T3D_fix3*, double, double)* @_ZNK8T3D_fix34callEdd to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.T3D_fix3*)* @_ZN8T3D_fix3D1Ev to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.T3D_fix3*)* @_ZN8T3D_fix3D0Ev to i32 (...)* (...)*) ], align 16, !dbg !2264
@_ZTV9Tinty_f2D = weak unnamed_addr global [5 x i32 (...)* (...)*]  [ i32 (...)* (...)*  null, i32 (...)* (...)*  bitcast (%struct.__si_class_type_info* @_ZTI9Tinty_f2D to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (double (%struct.Tinty_f2D*, double)* @_ZNK9Tinty_f2D4callEd to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.Tinty_f2D*)* @_ZN9Tinty_f2DD1Ev to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.Tinty_f2D*)* @_ZN9Tinty_f2DD0Ev to i32 (...)* (...)*) ], align 16, !dbg !1873
@_ZTV10Tintxy_f3D = weak unnamed_addr global [5 x i32 (...)* (...)*]  [ i32 (...)* (...)*  null, i32 (...)* (...)*  bitcast (%struct.__si_class_type_info* @_ZTI10Tintxy_f3D to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (double (%struct.Tintxy_f3D*, double)* @_ZNK10Tintxy_f3D4callEd to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.Tintxy_f3D*)* @_ZN10Tintxy_f3DD1Ev to i32 (...)* (...)*), i32 (...)* (...)*  bitcast (void (%struct.Tintxy_f3D*)* @_ZN10Tintxy_f3DD0Ev to i32 (...)* (...)*) ], align 16, !dbg !2012
@_ZTI11T1DFunction = weak unnamed_addr global %struct.__class_type_info  <{ %struct.__EDG_type_info  <{ i32 (...)* (...)*  bitcast (i8* getelementptr(i8, i8* bitcast([0 x i32 (...)* (...)*]* @_ZTVN10__cxxabiv117__class_type_infoE to i8*), i32 16) to i32 (...)* (...)*), i8*  getelementptr(i8, i8* bitcast([14 x i8]* @_ZTS11T1DFunction to i8*), i32 0) }> }>, align 16, !dbg !2837
@_ZTI11T2DFunction = weak unnamed_addr global %struct.__class_type_info  <{ %struct.__EDG_type_info  <{ i32 (...)* (...)*  bitcast (i8* getelementptr(i8, i8* bitcast([0 x i32 (...)* (...)*]* @_ZTVN10__cxxabiv117__class_type_infoE to i8*), i32 16) to i32 (...)* (...)*), i8*  getelementptr(i8, i8* bitcast([14 x i8]* @_ZTS11T2DFunction to i8*), i32 0) }> }>, align 16, !dbg !2839
@_ZTI8T2D_fix1 = weak unnamed_addr global %struct.__si_class_type_info  <{ %struct.__class_type_info  <{ %struct.__EDG_type_info  <{ i32 (...)* (...)*  bitcast (i8* getelementptr(i8, i8* bitcast([0 x i32 (...)* (...)*]* @_ZTVN10__cxxabiv120__si_class_type_infoE to i8*), i32 16) to i32 (...)* (...)*), i8*  getelementptr(i8, i8* bitcast([10 x i8]* @_ZTS8T2D_fix1 to i8*), i32 0) }> }>, %struct.__class_type_info*  @_ZTI11T1DFunction }>, align 16, !dbg !2841
@_ZTI8T3D_fix3 = weak unnamed_addr global %struct.__si_class_type_info  <{ %struct.__class_type_info  <{ %struct.__EDG_type_info  <{ i32 (...)* (...)*  bitcast (i8* getelementptr(i8, i8* bitcast([0 x i32 (...)* (...)*]* @_ZTVN10__cxxabiv120__si_class_type_infoE to i8*), i32 16) to i32 (...)* (...)*), i8*  getelementptr(i8, i8* bitcast([10 x i8]* @_ZTS8T3D_fix3 to i8*), i32 0) }> }>, %struct.__class_type_info*  @_ZTI11T2DFunction }>, align 16, !dbg !2843
@_ZTI9Tinty_f2D = weak unnamed_addr global %struct.__si_class_type_info  <{ %struct.__class_type_info  <{ %struct.__EDG_type_info  <{ i32 (...)* (...)*  bitcast (i8* getelementptr(i8, i8* bitcast([0 x i32 (...)* (...)*]* @_ZTVN10__cxxabiv120__si_class_type_infoE to i8*), i32 16) to i32 (...)* (...)*), i8*  getelementptr(i8, i8* bitcast([11 x i8]* @_ZTS9Tinty_f2D to i8*), i32 0) }> }>, %struct.__class_type_info*  @_ZTI11T1DFunction }>, align 16, !dbg !2845
@_ZTI10Tintxy_f3D = weak unnamed_addr global %struct.__si_class_type_info  <{ %struct.__class_type_info  <{ %struct.__EDG_type_info  <{ i32 (...)* (...)*  bitcast (i8* getelementptr(i8, i8* bitcast([0 x i32 (...)* (...)*]* @_ZTVN10__cxxabiv120__si_class_type_infoE to i8*), i32 16) to i32 (...)* (...)*), i8*  getelementptr(i8, i8* bitcast([13 x i8]* @_ZTS10Tintxy_f3D to i8*), i32 0) }> }>, %struct.__class_type_info*  @_ZTI11T1DFunction }>, align 16, !dbg !2847
@_ZTS11T1DFunction = weak unnamed_addr global [14 x i8]  c"11T1DFunction\00", align 8, !dbg !2855
@_ZTS11T2DFunction = weak unnamed_addr global [14 x i8]  c"11T2DFunction\00", align 8, !dbg !2857
@_ZTS8T2D_fix1 = weak unnamed_addr global [10 x i8]  c"8T2D_fix1\00", align 8, !dbg !2864
@_ZTS8T3D_fix3 = weak unnamed_addr global [10 x i8]  c"8T3D_fix3\00", align 8, !dbg !2866
@_ZTS9Tinty_f2D = weak unnamed_addr global [11 x i8]  c"9Tinty_f2D\00", align 8, !dbg !2871
@_ZTS10Tintxy_f3D = weak unnamed_addr global [13 x i8]  c"10Tintxy_f3D\00", align 8, !dbg !2874
@__I___25_backgroundfield_quadr_cpp_6e0a9365 = global i32 0, align 4, !dbg !2828
@_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE = internal global %struct._ZNSt8ios_base4InitE zeroinitializer , align 1, !dbg !2833
@__dso_handle = external global i8*, align 8
@_ZSt4cerr = external global %struct._ZSo, align 8
@_ZTVN10__cxxabiv120__si_class_type_infoE = external global [0 x i32 (...)* (...)*], align 8
@_ZTVN10__cxxabiv117__class_type_infoE = external global [0 x i32 (...)* (...)*], align 8
@.S08003 = internal constant [21 x i8] [i8 42,i8 42,i8 42,i8 32,i8 69,i8 114,i8 114,i8 111,i8 114,i8 32,i8 105,i8 110,i8 32,i8 114,i8 97,i8 116,i8 105,i8 110,i8 116,i8 10,i8 0], align 1
@llvm.global_ctors = appending global [1 x { i32, void ()*, i8* }][{ i32, void ()*, i8* } { i32 65535, void ()* @__sti___25_backgroundfield_quadr_cpp_6e0a9365, i8* null }]
attributes #0 = { "frame-pointer"="all" }

declare void @__cxa_pure_virtual() #0
declare signext i32 @__cxa_atexit(void (i8*)*, i8*, i8*) #0
declare void @_ZNSt8ios_base4InitD1Ev(%struct._ZNSt8ios_base4InitE*) #0
declare void @_ZNSt8ios_base4InitC1Ev(%struct._ZNSt8ios_base4InitE*) #0
declare void @_ZdlPvm(i8*, i64) #0
declare double @llvm.maxnum.f64(double, double)
declare void @exit(i32 signext) #0
declare %struct._ZSo* @_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l(%struct._ZSo*, i8*, i64) #0
declare i64 @strlen(i8*) #0
declare void @_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate(%struct._ZSt9basic_iosIcSt11char_traitsIcEE*, i32 signext) #0
declare double @llvm.fabs.f64(double)
declare double @llvm.fma.f64(double, double, double)
declare void @llvm.dbg.declare(metadata, metadata, metadata)
declare i32 @__gxx_personality_v0(...)

; Named metadata
!llvm.module.flags = !{ !1, !2 }
!llvm.dbg.cu = !{ !10 }

; Metadata
!1 = !{ i32 2, !"Dwarf Version", i32 2 }
!2 = !{ i32 2, !"Debug Info Version", i32 3 }
!3 = !DIFile(filename: "backgroundfield/quadr.cpp", directory: "/home/talgat/vlasiator")
; !4 = !DIFile(tag: DW_TAG_file_type, pair: !3)
!4 = !{ i32 41, !3 }
!5 = !{  }
!6 = !{  }
!7 = !{ !1359, !1408, !1453, !1537, !1616, !1761, !1901, !2039, !2047, !2057, !2067, !2075, !2088, !2098, !2108, !2118, !2126, !2138, !2160, !2178, !2190, !2206, !2224, !2244, !2266, !2284, !2298, !2314, !2332, !2352, !2376, !2394, !2515, !2532, !2550, !2570, !2598, !2616, !2666, !2683, !2701, !2723, !2737, !2747, !2770, !2789, !2822 }
!8 = !{ !1524, !1871, !1873, !2012, !2096, !2158, !2264, !2828, !2833, !2835, !2837, !2839, !2841, !2843, !2845, !2847, !2850, !2855, !2857, !2859, !2864, !2866, !2871, !2874 }
!9 = !{  }
!10 = distinct !DICompileUnit(file: !3, language: DW_LANG_C_plus_plus, producer: " NVC++ 21.2-0", enums: !5, retainedTypes: !6, globals: !8, emissionKind: FullDebug, imports: !9)
!11 = !DINamespace(scope: !10, name: "std")
!12 = !DINamespace(scope: !11, name: "__cxx11")
!13 = !DINamespace(scope: !11, name: "__exception_ptr")
!14 = !{ !18 }
!15 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !13, name: "exception_ptr", line: 1, size: 64, align: 64, elements: !14, runtimeLang: DW_LANG_C_plus_plus)
!16 = !DIBasicType(tag: DW_TAG_unspecified_type, name: "void")
!17 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !16)
!18 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !15, name: "_M_exception_object", line: 1, size: 64, align: 64, baseType: !17)
!19 = !DINamespace(scope: !11, name: "__swappable_details")
!20 = !{  }
!21 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !19, name: "__do_is_swappable_impl", line: 1, size: 8, align: 8, elements: !20, runtimeLang: DW_LANG_C_plus_plus)
!22 = !{  }
!23 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !19, name: "__do_is_nothrow_swappable_impl", line: 1, size: 8, align: 8, elements: !22, runtimeLang: DW_LANG_C_plus_plus)
!24 = !DINamespace(scope: !11, name: "__swappable_with_details")
!25 = !{  }
!26 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !24, name: "__do_is_swappable_with_impl", line: 1, size: 8, align: 8, elements: !25, runtimeLang: DW_LANG_C_plus_plus)
!27 = !{  }
!28 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !24, name: "__do_is_nothrow_swappable_with_impl", line: 1, size: 8, align: 8, elements: !27, runtimeLang: DW_LANG_C_plus_plus)
!29 = !DINamespace(scope: !11, name: "__debug")
!30 = !DIBasicType(tag: DW_TAG_base_type, name: "unsigned long", size: 64, align: 64, encoding: DW_ATE_unsigned)
!31 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "size_t", line: 1, size: 64, align: 64, baseType: !30)
!32 = !DIBasicType(tag: DW_TAG_base_type, name: "long", size: 64, align: 64, encoding: DW_ATE_signed)
!33 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "ptrdiff_t", line: 1, size: 64, align: 64, baseType: !32)
!34 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "nullptr_t", line: 1, size: 64, align: 64, baseType: !17)
!35 = !{  }
!36 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__true_type", line: 1, size: 8, align: 8, elements: !35, runtimeLang: DW_LANG_C_plus_plus)
!37 = !{  }
!38 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__false_type", line: 1, size: 8, align: 8, elements: !37, runtimeLang: DW_LANG_C_plus_plus)
!39 = !{ !46, !47, !55 }
!40 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE", size: 256, align: 64, elements: !39)
!41 = !{ !45 }
!42 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_Alloc_hiderE", size: 64, align: 64, elements: !41)
!43 = !DIBasicType(tag: DW_TAG_base_type, name: "signed char", size: 8, align: 8, encoding: DW_ATE_signed_char)
!44 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !43)
!45 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !42, name: "_M_p", size: 64, align: 64, baseType: !44)
!46 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !40, name: "_M_dataplus", size: 64, align: 64, baseType: !42)
!47 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !40, name: "_M_string_length", size: 64, align: 64, offset: 64, baseType: !30)
!48 = !{ !53, !54 }
!49 = !DICompositeType(tag: DW_TAG_union_type, file: !3, name: "__C5", size: 128, align: 64, elements: !48)
!50 = !DISubrange(count: 16)
!51 = !{ !50 }
!52 = !DICompositeType(tag: DW_TAG_array_type, size: 128, align: 8, baseType: !43, elements: !51)
!53 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !49, name: "_M_local_buf", size: 128, align: 8, baseType: !52)
!54 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !49, name: "_M_allocated_capacity", size: 64, align: 64, baseType: !30)
!55 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !40, size: 128, align: 64, offset: 128, baseType: !49)
!56 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "string", line: 1, size: 256, align: 64, baseType: !40)
!57 = !{ !64, !65, !73 }
!58 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt7__cxx1112basic_stringIwSt11char_traitsIwESaIwEEE", size: 256, align: 64, elements: !57)
!59 = !{ !63 }
!60 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt7__cxx1112basic_stringIwSt11char_traitsIwESaIwEE12_Alloc_hiderE", size: 64, align: 64, elements: !59)
!61 = !DIBasicType(tag: DW_TAG_base_type, name: "int", size: 32, align: 32, encoding: DW_ATE_signed)
!62 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !61)
!63 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !60, name: "_M_p", size: 64, align: 64, baseType: !62)
!64 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !58, name: "_M_dataplus", size: 64, align: 64, baseType: !60)
!65 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !58, name: "_M_string_length", size: 64, align: 64, offset: 64, baseType: !30)
!66 = !{ !71, !72 }
!67 = !DICompositeType(tag: DW_TAG_union_type, file: !3, name: "__C6", size: 128, align: 64, elements: !66)
!68 = !DISubrange(count: 4)
!69 = !{ !68 }
!70 = !DICompositeType(tag: DW_TAG_array_type, size: 128, align: 32, baseType: !61, elements: !69)
!71 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !67, name: "_M_local_buf", size: 128, align: 32, baseType: !70)
!72 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !67, name: "_M_allocated_capacity", size: 64, align: 64, baseType: !30)
!73 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !58, size: 128, align: 64, offset: 128, baseType: !67)
!74 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "wstring", line: 1, size: 256, align: 64, baseType: !58)
!75 = !{ !82, !83, !91 }
!76 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt7__cxx1112basic_stringIDsSt11char_traitsIDsESaIDsEEE", size: 256, align: 64, elements: !75)
!77 = !{ !81 }
!78 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt7__cxx1112basic_stringIDsSt11char_traitsIDsESaIDsEE12_Alloc_hiderE", size: 64, align: 64, elements: !77)
!79 = !DIBasicType(tag: DW_TAG_base_type, name: "unsigned short", size: 16, align: 16, encoding: DW_ATE_unsigned)
!80 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !79)
!81 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !78, name: "_M_p", size: 64, align: 64, baseType: !80)
!82 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !76, name: "_M_dataplus", size: 64, align: 64, baseType: !78)
!83 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !76, name: "_M_string_length", size: 64, align: 64, offset: 64, baseType: !30)
!84 = !{ !89, !90 }
!85 = !DICompositeType(tag: DW_TAG_union_type, file: !3, name: "__C7", size: 128, align: 64, elements: !84)
!86 = !DISubrange(count: 8)
!87 = !{ !86 }
!88 = !DICompositeType(tag: DW_TAG_array_type, size: 128, align: 16, baseType: !79, elements: !87)
!89 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !85, name: "_M_local_buf", size: 128, align: 16, baseType: !88)
!90 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !85, name: "_M_allocated_capacity", size: 64, align: 64, baseType: !30)
!91 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !76, size: 128, align: 64, offset: 128, baseType: !85)
!92 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "u16string", line: 1, size: 256, align: 64, baseType: !76)
!93 = !{ !100, !101, !107 }
!94 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt7__cxx1112basic_stringIDiSt11char_traitsIDiESaIDiEEE", size: 256, align: 64, elements: !93)
!95 = !{ !99 }
!96 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt7__cxx1112basic_stringIDiSt11char_traitsIDiESaIDiEE12_Alloc_hiderE", size: 64, align: 64, elements: !95)
!97 = !DIBasicType(tag: DW_TAG_base_type, name: "unsigned", size: 32, align: 32, encoding: DW_ATE_unsigned)
!98 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !97)
!99 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !96, name: "_M_p", size: 64, align: 64, baseType: !98)
!100 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !94, name: "_M_dataplus", size: 64, align: 64, baseType: !96)
!101 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !94, name: "_M_string_length", size: 64, align: 64, offset: 64, baseType: !30)
!102 = !{ !105, !106 }
!103 = !DICompositeType(tag: DW_TAG_union_type, file: !3, name: "__C8", size: 128, align: 64, elements: !102)
!104 = !DICompositeType(tag: DW_TAG_array_type, size: 128, align: 32, baseType: !97, elements: !69)
!105 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !103, name: "_M_local_buf", size: 128, align: 32, baseType: !104)
!106 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !103, name: "_M_allocated_capacity", size: 64, align: 64, baseType: !30)
!107 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !94, size: 128, align: 64, offset: 128, baseType: !103)
!108 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "u32string", line: 1, size: 256, align: 64, baseType: !94)
!109 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "streamoff", line: 1, size: 64, align: 64, baseType: !32)
!110 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "streamsize", line: 1, size: 64, align: 64, baseType: !32)
!111 = !{  }
!112 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt4fposI11__mbstate_tE", align: 8, elements: !111)
!113 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "streampos", line: 1, align: 8, baseType: !112)
!114 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "wstreampos", line: 1, align: 8, baseType: !112)
!115 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "u16streampos", line: 1, align: 8, baseType: !112)
!116 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "u32streampos", line: 1, align: 8, baseType: !112)
!117 = !{ !211, !217, !218, !219, !231, !267, !272, !277 }
!118 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt9basic_iosIcSt11char_traitsIcEE", size: 2112, align: 64, elements: !117)
!119 = !{ !127, !128, !129, !153, !163, !164, !181, !186, !188, !189, !191, !210 }
!120 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt8ios_base", size: 1728, align: 64, elements: !119)
!121 = !{ !61 }
!122 = !DISubroutineType(types: !121)
!123 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !122)
!124 = !{ !123 }
!125 = !DISubroutineType(types: !124)
!126 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !125)
!127 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "__vptr", size: 64, align: 64, baseType: !126)
!128 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_precision", size: 64, align: 64, offset: 64, baseType: !32)
!129 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_width", size: 64, align: 64, offset: 128, baseType: !32)
!130 = !DIEnumerator(name: "_ZSt19_S_ios_fmtflags_min", value: -2147483648)
!131 = !DIEnumerator(name: "_ZSt19_S_ios_fmtflags_max", value: 2147483647)
!132 = !DIEnumerator(name: "_ZSt19_S_ios_fmtflags_end", value: 65536)
!133 = !DIEnumerator(name: "_ZSt13_S_floatfield", value: 260)
!134 = !DIEnumerator(name: "_ZSt12_S_basefield", value: 74)
!135 = !DIEnumerator(name: "_ZSt14_S_adjustfield", value: 176)
!136 = !DIEnumerator(name: "_ZSt12_S_uppercase", value: 16384)
!137 = !DIEnumerator(name: "_ZSt10_S_unitbuf", value: 8192)
!138 = !DIEnumerator(name: "_ZSt9_S_skipws", value: 4096)
!139 = !DIEnumerator(name: "_ZSt10_S_showpos", value: 2048)
!140 = !DIEnumerator(name: "_ZSt12_S_showpoint", value: 1024)
!141 = !DIEnumerator(name: "_ZSt11_S_showbase", value: 512)
!142 = !DIEnumerator(name: "_ZSt13_S_scientific", value: 256)
!143 = !DIEnumerator(name: "_ZSt8_S_right", value: 128)
!144 = !DIEnumerator(name: "_ZSt6_S_oct", value: 64)
!145 = !DIEnumerator(name: "_ZSt7_S_left", value: 32)
!146 = !DIEnumerator(name: "_ZSt11_S_internal", value: 16)
!147 = !DIEnumerator(name: "_ZSt6_S_hex", value: 8)
!148 = !DIEnumerator(name: "_ZSt8_S_fixed", value: 4)
!149 = !DIEnumerator(name: "_ZSt6_S_dec", value: 2)
!150 = !DIEnumerator(name: "_ZSt12_S_boolalpha", value: 1)
!151 = !{ !150, !149, !148, !147, !146, !145, !144, !143, !142, !141, !140, !139, !138, !137, !136, !135, !134, !133, !132, !131, !130 }
!152 = !DICompositeType(tag: DW_TAG_enumeration_type, file: !3, name: "_ZSt13_Ios_Fmtflags", size: 32, align: 32, elements: !151)
!153 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_flags", size: 32, align: 32, offset: 192, baseType: !152)
!154 = !DIEnumerator(name: "_ZSt18_S_ios_iostate_min", value: -2147483648)
!155 = !DIEnumerator(name: "_ZSt18_S_ios_iostate_max", value: 2147483647)
!156 = !DIEnumerator(name: "_ZSt18_S_ios_iostate_end", value: 65536)
!157 = !DIEnumerator(name: "_ZSt10_S_failbit", value: 4)
!158 = !DIEnumerator(name: "_ZSt9_S_eofbit", value: 2)
!159 = !DIEnumerator(name: "_ZSt9_S_badbit", value: 1)
!160 = !DIEnumerator(name: "_ZSt10_S_goodbit", value: 0)
!161 = !{ !160, !159, !158, !157, !156, !155, !154 }
!162 = !DICompositeType(tag: DW_TAG_enumeration_type, file: !3, name: "_ZSt12_Ios_Iostate", size: 32, align: 32, elements: !161)
!163 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_exception", size: 32, align: 32, offset: 224, baseType: !162)
!164 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_streambuf_state", size: 32, align: 32, offset: 256, baseType: !162)
!165 = !{ !168, !178, !179, !180 }
!166 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt8ios_base14_Callback_listE", size: 192, align: 64, elements: !165)
!167 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !166)
!168 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !166, name: "_M_next", size: 64, align: 64, baseType: !167)
!169 = !DIEnumerator(name: "_ZNSt8ios_base13copyfmt_eventE", value: 2)
!170 = !DIEnumerator(name: "_ZNSt8ios_base11imbue_eventE", value: 1)
!171 = !DIEnumerator(name: "_ZNSt8ios_base11erase_eventE", value: 0)
!172 = !{ !171, !170, !169 }
!173 = !DICompositeType(tag: DW_TAG_enumeration_type, file: !3, name: "_ZNSt8ios_base5eventE", size: 32, align: 32, elements: !172)
!174 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !120)
!175 = !{ null, !173, !174, !61 }
!176 = !DISubroutineType(types: !175)
!177 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !176)
!178 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !166, name: "_M_fn", size: 64, align: 64, offset: 64, baseType: !177)
!179 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !166, name: "_M_index", size: 32, align: 32, offset: 128, baseType: !61)
!180 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !166, name: "_M_refcount", size: 32, align: 32, offset: 160, baseType: !61)
!181 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_callbacks", size: 64, align: 64, offset: 320, baseType: !167)
!182 = !{ !184, !185 }
!183 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt8ios_base6_WordsE", size: 128, align: 64, elements: !182)
!184 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !183, name: "_M_pword", size: 64, align: 64, baseType: !17)
!185 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !183, name: "_M_iword", size: 64, align: 64, offset: 64, baseType: !32)
!186 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_word_zero", size: 128, align: 64, offset: 384, baseType: !183)
!187 = !DICompositeType(tag: DW_TAG_array_type, size: 1024, align: 64, baseType: !183, elements: !87)
!188 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_local_word", size: 1024, align: 64, offset: 512, baseType: !187)
!189 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_word_size", size: 32, align: 32, offset: 1536, baseType: !61)
!190 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !183)
!191 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_word", size: 64, align: 64, offset: 1600, baseType: !190)
!192 = !{ !209 }
!193 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt6locale", size: 64, align: 64, elements: !192)
!194 = !{ !196, !203, !204, !205, !207 }
!195 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt6locale5_ImplE", size: 320, align: 64, elements: !194)
!196 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !195, name: "_M_refcount", size: 32, align: 32, baseType: !61)
!197 = !{ !199, !200 }
!198 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt6locale5facetE", size: 128, align: 64, elements: !197)
!199 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !198, name: "__vptr", size: 64, align: 64, baseType: !126)
!200 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !198, name: "_M_refcount", size: 32, align: 32, offset: 64, baseType: !61)
!201 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !198)
!202 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !201)
!203 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !195, name: "_M_facets", size: 64, align: 64, offset: 64, baseType: !202)
!204 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !195, name: "_M_facets_size", size: 64, align: 64, offset: 128, baseType: !30)
!205 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !195, name: "_M_caches", size: 64, align: 64, offset: 192, baseType: !202)
!206 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !44)
!207 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !195, name: "_M_names", size: 64, align: 64, offset: 256, baseType: !206)
!208 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !195)
!209 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !193, name: "_M_impl", size: 64, align: 64, baseType: !208)
!210 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !120, name: "_M_ios_locale", size: 64, align: 64, offset: 1664, baseType: !193)
!211 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !118, name: "__b_St8ios_base", size: 1728, align: 64, baseType: !120)
!212 = !{ !214, !215 }
!213 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSo", size: 2176, align: 64, elements: !212)
!214 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !213, name: "__vptr", size: 64, align: 64, baseType: !126)
!215 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !213, name: "__v_St9basic_iosIcSt11char_traitsIcEE", size: 2112, align: 64, offset: 64, baseType: !118)
!216 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !213)
!217 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !118, name: "_M_tie", size: 64, align: 64, offset: 1728, baseType: !216)
!218 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !118, name: "_M_fill", size: 8, align: 8, offset: 1792, baseType: !43)
!219 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !118, name: "_M_fill_init", size: 8, align: 8, offset: 1800, baseType: !43)
!220 = !{ !222, !223, !224, !225, !226, !227, !228, !229 }
!221 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt15basic_streambufIcSt11char_traitsIcEE", size: 512, align: 64, elements: !220)
!222 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !221, name: "__vptr", size: 64, align: 64, baseType: !126)
!223 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !221, name: "_M_in_beg", size: 64, align: 64, offset: 64, baseType: !44)
!224 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !221, name: "_M_in_cur", size: 64, align: 64, offset: 128, baseType: !44)
!225 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !221, name: "_M_in_end", size: 64, align: 64, offset: 192, baseType: !44)
!226 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !221, name: "_M_out_beg", size: 64, align: 64, offset: 256, baseType: !44)
!227 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !221, name: "_M_out_cur", size: 64, align: 64, offset: 320, baseType: !44)
!228 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !221, name: "_M_out_end", size: 64, align: 64, offset: 384, baseType: !44)
!229 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !221, name: "_M_buf_locale", size: 64, align: 64, offset: 448, baseType: !193)
!230 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !221)
!231 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !118, name: "_M_streambuf", size: 64, align: 64, offset: 1856, baseType: !230)
!232 = !{ !238, !254, !255, !256, !257, !258, !259, !263, !264, !265 }
!233 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt5ctypeIcE", size: 4608, align: 64, elements: !232)
!234 = !{ !236, !237 }
!235 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "__SO__NSt6locale5facetE", size: 96, align: 64, elements: !234)
!236 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !235, name: "__vptr", size: 64, align: 64, baseType: !126)
!237 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !235, name: "_M_refcount", size: 32, align: 32, offset: 64, baseType: !61)
!238 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "__b_NSt6locale5facetE", size: 96, align: 64, baseType: !235)
!239 = !{ !247, !248, !249, !250, !252 }
!240 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "__locale_struct", size: 1856, align: 64, elements: !239)
!241 = !DISubrange(count: 13)
!242 = !{  }
!243 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "__locale_data", align: 8, elements: !242)
!244 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !243)
!245 = !{ !241 }
!246 = !DICompositeType(tag: DW_TAG_array_type, size: 832, align: 64, baseType: !244, elements: !245)
!247 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !240, name: "__locales", size: 832, align: 64, baseType: !246)
!248 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !240, name: "__ctype_b", size: 64, align: 64, offset: 832, baseType: !80)
!249 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !240, name: "__ctype_tolower", size: 64, align: 64, offset: 896, baseType: !62)
!250 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !240, name: "__ctype_toupper", size: 64, align: 64, offset: 960, baseType: !62)
!251 = !DICompositeType(tag: DW_TAG_array_type, size: 832, align: 64, baseType: !44, elements: !245)
!252 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !240, name: "__names", size: 832, align: 64, offset: 1024, baseType: !251)
!253 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !240)
!254 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_c_locale_ctype", size: 64, align: 64, offset: 128, baseType: !253)
!255 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_del", size: 8, align: 8, offset: 192, baseType: !43)
!256 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_toupper", size: 64, align: 64, offset: 256, baseType: !62)
!257 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_tolower", size: 64, align: 64, offset: 320, baseType: !62)
!258 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_table", size: 64, align: 64, offset: 384, baseType: !80)
!259 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_widen_ok", size: 8, align: 8, offset: 448, baseType: !43)
!260 = !DISubrange(count: 256)
!261 = !{ !260 }
!262 = !DICompositeType(tag: DW_TAG_array_type, size: 2048, align: 8, baseType: !43, elements: !261)
!263 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_widen", size: 2048, align: 8, offset: 456, baseType: !262)
!264 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_narrow", size: 2048, align: 8, offset: 2504, baseType: !262)
!265 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !233, name: "_M_narrow_ok", size: 8, align: 8, offset: 4552, baseType: !43)
!266 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !233)
!267 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !118, name: "_M_ctype", size: 64, align: 64, offset: 1920, baseType: !266)
!268 = !{ !270 }
!269 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt7num_putIcSt19ostreambuf_iteratorIcSt11char_traitsIcEEE", size: 128, align: 64, elements: !268)
!270 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !269, name: "__b_NSt6locale5facetE", size: 96, align: 64, baseType: !235)
!271 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !269)
!272 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !118, name: "_M_num_put", size: 64, align: 64, offset: 1984, baseType: !271)
!273 = !{ !275 }
!274 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt7num_getIcSt19istreambuf_iteratorIcSt11char_traitsIcEEE", size: 128, align: 64, elements: !273)
!275 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !274, name: "__b_NSt6locale5facetE", size: 96, align: 64, baseType: !235)
!276 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !274)
!277 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !118, name: "_M_num_get", size: 64, align: 64, offset: 2048, baseType: !276)
!278 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "ios", line: 1, size: 2112, align: 64, baseType: !118)
!279 = !{ !281, !282, !283 }
!280 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSi", size: 2240, align: 64, elements: !279)
!281 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !280, name: "__vptr", size: 64, align: 64, baseType: !126)
!282 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !280, name: "_M_gcount", size: 64, align: 64, offset: 64, baseType: !32)
!283 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !280, name: "__v_St9basic_iosIcSt11char_traitsIcEE", size: 2112, align: 64, offset: 128, baseType: !118)
!284 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "istream", line: 1, size: 2240, align: 64, baseType: !280)
!285 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "ostream", line: 1, size: 2176, align: 64, baseType: !213)
!286 = !{ !288, !289, !343 }
!287 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt13basic_istreamIwSt11char_traitsIwEE", size: 2240, align: 64, elements: !286)
!288 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !287, name: "__vptr", size: 64, align: 64, baseType: !126)
!289 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !287, name: "_M_gcount", size: 64, align: 64, offset: 64, baseType: !32)
!290 = !{ !292, !298, !299, !300, !312, !332, !337, !342 }
!291 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt9basic_iosIwSt11char_traitsIwEE", size: 2112, align: 64, elements: !290)
!292 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !291, name: "__b_St8ios_base", size: 1728, align: 64, baseType: !120)
!293 = !{ !295, !296 }
!294 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt13basic_ostreamIwSt11char_traitsIwEE", size: 2176, align: 64, elements: !293)
!295 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !294, name: "__vptr", size: 64, align: 64, baseType: !126)
!296 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !294, name: "__v_St9basic_iosIwSt11char_traitsIwEE", size: 2112, align: 64, offset: 64, baseType: !291)
!297 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !294)
!298 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !291, name: "_M_tie", size: 64, align: 64, offset: 1728, baseType: !297)
!299 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !291, name: "_M_fill", size: 32, align: 32, offset: 1792, baseType: !61)
!300 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !291, name: "_M_fill_init", size: 8, align: 8, offset: 1824, baseType: !43)
!301 = !{ !303, !304, !305, !306, !307, !308, !309, !310 }
!302 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt15basic_streambufIwSt11char_traitsIwEE", size: 512, align: 64, elements: !301)
!303 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !302, name: "__vptr", size: 64, align: 64, baseType: !126)
!304 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !302, name: "_M_in_beg", size: 64, align: 64, offset: 64, baseType: !62)
!305 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !302, name: "_M_in_cur", size: 64, align: 64, offset: 128, baseType: !62)
!306 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !302, name: "_M_in_end", size: 64, align: 64, offset: 192, baseType: !62)
!307 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !302, name: "_M_out_beg", size: 64, align: 64, offset: 256, baseType: !62)
!308 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !302, name: "_M_out_cur", size: 64, align: 64, offset: 320, baseType: !62)
!309 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !302, name: "_M_out_end", size: 64, align: 64, offset: 384, baseType: !62)
!310 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !302, name: "_M_buf_locale", size: 64, align: 64, offset: 448, baseType: !193)
!311 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !302)
!312 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !291, name: "_M_streambuf", size: 64, align: 64, offset: 1856, baseType: !311)
!313 = !{ !318, !319, !320, !324, !326, !328, !330 }
!314 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt5ctypeIwE", size: 10752, align: 64, elements: !313)
!315 = !{ !317 }
!316 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "__SO__St21__ctype_abstract_baseIwE", size: 96, align: 64, elements: !315)
!317 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !316, name: "__b_NSt6locale5facetE", size: 96, align: 64, baseType: !235)
!318 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !314, name: "__b_St21__ctype_abstract_baseIwE", size: 96, align: 64, baseType: !316)
!319 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !314, name: "_M_c_locale_ctype", size: 64, align: 64, offset: 128, baseType: !253)
!320 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !314, name: "_M_narrow_ok", size: 8, align: 8, offset: 192, baseType: !43)
!321 = !DISubrange(count: 128)
!322 = !{ !321 }
!323 = !DICompositeType(tag: DW_TAG_array_type, size: 1024, align: 8, baseType: !43, elements: !322)
!324 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !314, name: "_M_narrow", size: 1024, align: 8, offset: 200, baseType: !323)
!325 = !DICompositeType(tag: DW_TAG_array_type, size: 8192, align: 32, baseType: !97, elements: !261)
!326 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !314, name: "_M_widen", size: 8192, align: 32, offset: 1248, baseType: !325)
!327 = !DICompositeType(tag: DW_TAG_array_type, size: 256, align: 16, baseType: !79, elements: !51)
!328 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !314, name: "_M_bit", size: 256, align: 16, offset: 9440, baseType: !327)
!329 = !DICompositeType(tag: DW_TAG_array_type, size: 1024, align: 64, baseType: !30, elements: !51)
!330 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !314, name: "_M_wmask", size: 1024, align: 64, offset: 9728, baseType: !329)
!331 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !314)
!332 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !291, name: "_M_ctype", size: 64, align: 64, offset: 1920, baseType: !331)
!333 = !{ !335 }
!334 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt7num_putIwSt19ostreambuf_iteratorIwSt11char_traitsIwEEE", size: 128, align: 64, elements: !333)
!335 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !334, name: "__b_NSt6locale5facetE", size: 96, align: 64, baseType: !235)
!336 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !334)
!337 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !291, name: "_M_num_put", size: 64, align: 64, offset: 1984, baseType: !336)
!338 = !{ !340 }
!339 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt7num_getIwSt19istreambuf_iteratorIwSt11char_traitsIwEEE", size: 128, align: 64, elements: !338)
!340 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !339, name: "__b_NSt6locale5facetE", size: 96, align: 64, baseType: !235)
!341 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !339)
!342 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !291, name: "_M_num_get", size: 64, align: 64, offset: 2048, baseType: !341)
!343 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !287, name: "__v_St9basic_iosIwSt11char_traitsIwEE", size: 2112, align: 64, offset: 128, baseType: !291)
!344 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "wistream", line: 1, size: 2240, align: 64, baseType: !287)
!345 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "wostream", line: 1, size: 2176, align: 64, baseType: !294)
!346 = !{  }
!347 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "exception", line: 1, size: 64, align: 64, elements: !346, runtimeLang: DW_LANG_C_plus_plus)
!348 = !{ !350 }
!349 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "bad_exception", line: 1, size: 64, align: 64, elements: !348, runtimeLang: DW_LANG_C_plus_plus)
!350 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !349, name: "exception", line: 1, size: 64, align: 64, baseType: !347)
!351 = !{ null }
!352 = !DISubroutineType(types: !351)
!353 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !352)
!354 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "terminate_handler", line: 1, size: 64, align: 64, baseType: !353)
!355 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "unexpected_handler", line: 1, size: 64, align: 64, baseType: !353)
!356 = !{ !358 }
!357 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "type_info", line: 1, size: 128, align: 64, elements: !356, runtimeLang: DW_LANG_C_plus_plus)
!358 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !357, name: "__name", line: 1, size: 64, align: 64, offset: 64, baseType: !44)
!359 = !{ !361 }
!360 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "bad_cast", line: 1, size: 64, align: 64, elements: !359, runtimeLang: DW_LANG_C_plus_plus)
!361 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !360, name: "exception", line: 1, size: 64, align: 64, baseType: !347)
!362 = !{ !364 }
!363 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "bad_typeid", line: 1, size: 64, align: 64, elements: !362, runtimeLang: DW_LANG_C_plus_plus)
!364 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !363, name: "exception", line: 1, size: 64, align: 64, baseType: !347)
!365 = !{ !367 }
!366 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "bad_alloc", line: 1, size: 64, align: 64, elements: !365, runtimeLang: DW_LANG_C_plus_plus)
!367 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !366, name: "exception", line: 1, size: 64, align: 64, baseType: !347)
!368 = !{ !370 }
!369 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "bad_array_new_length", line: 1, size: 64, align: 64, elements: !368, runtimeLang: DW_LANG_C_plus_plus)
!370 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !369, name: "bad_alloc", line: 1, size: 64, align: 64, baseType: !366)
!371 = !{  }
!372 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "nothrow_t", line: 1, size: 8, align: 8, elements: !371, runtimeLang: DW_LANG_C_plus_plus)
!373 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "new_handler", line: 1, size: 64, align: 64, baseType: !353)
!374 = !{ !379 }
!375 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt17integral_constantIbLb1EE", size: 8, align: 8, elements: !374)
!376 = !DISubrange(count: 0)
!377 = !{ !376 }
!378 = !DICompositeType(tag: DW_TAG_array_type, size: 8, align: 8, baseType: !43, elements: !377)
!379 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !375, size: 8, align: 8, baseType: !378)
!380 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "true_type", line: 1, size: 8, align: 8, baseType: !375)
!381 = !{ !383 }
!382 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZSt17integral_constantIbLb0EE", size: 8, align: 8, elements: !381)
!383 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !382, size: 8, align: 8, baseType: !378)
!384 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "false_type", line: 1, size: 8, align: 8, baseType: !382)
!385 = !{  }
!386 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__failure_type", line: 1, size: 8, align: 8, elements: !385, runtimeLang: DW_LANG_C_plus_plus)
!387 = !{  }
!388 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__do_is_destructible_impl", line: 1, size: 8, align: 8, elements: !387, runtimeLang: DW_LANG_C_plus_plus)
!389 = !{  }
!390 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__do_is_nt_destructible_impl", line: 1, size: 8, align: 8, elements: !389, runtimeLang: DW_LANG_C_plus_plus)
!391 = !{  }
!392 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__do_is_implicitly_default_constructible_impl", line: 1, size: 8, align: 8, elements: !391, runtimeLang: DW_LANG_C_plus_plus)
!393 = !{  }
!394 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "__make_unsigned_selector_base", line: 1, size: 8, align: 8, elements: !393, runtimeLang: DW_LANG_C_plus_plus)
!395 = !{  }
!396 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__do_common_type_impl", line: 1, size: 8, align: 8, elements: !395, runtimeLang: DW_LANG_C_plus_plus)
!397 = !{  }
!398 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__do_member_type_wrapper", line: 1, size: 8, align: 8, elements: !397, runtimeLang: DW_LANG_C_plus_plus)
!399 = !{  }
!400 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__invoke_memfun_ref", line: 1, size: 8, align: 8, elements: !399, runtimeLang: DW_LANG_C_plus_plus)
!401 = !{  }
!402 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__invoke_memfun_deref", line: 1, size: 8, align: 8, elements: !401, runtimeLang: DW_LANG_C_plus_plus)
!403 = !{  }
!404 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__invoke_memobj_ref", line: 1, size: 8, align: 8, elements: !403, runtimeLang: DW_LANG_C_plus_plus)
!405 = !{  }
!406 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__invoke_memobj_deref", line: 1, size: 8, align: 8, elements: !405, runtimeLang: DW_LANG_C_plus_plus)
!407 = !{  }
!408 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__invoke_other", line: 1, size: 8, align: 8, elements: !407, runtimeLang: DW_LANG_C_plus_plus)
!409 = !{  }
!410 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__result_of_memfun_ref_impl", line: 1, size: 8, align: 8, elements: !409, runtimeLang: DW_LANG_C_plus_plus)
!411 = !{  }
!412 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__result_of_memfun_deref_impl", line: 1, size: 8, align: 8, elements: !411, runtimeLang: DW_LANG_C_plus_plus)
!413 = !{  }
!414 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__result_of_memobj_ref_impl", line: 1, size: 8, align: 8, elements: !413, runtimeLang: DW_LANG_C_plus_plus)
!415 = !{  }
!416 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__result_of_memobj_deref_impl", line: 1, size: 8, align: 8, elements: !415, runtimeLang: DW_LANG_C_plus_plus)
!417 = !{  }
!418 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__result_of_other_impl", line: 1, size: 8, align: 8, elements: !417, runtimeLang: DW_LANG_C_plus_plus)
!419 = !{  }
!420 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__nonesuch", line: 1, size: 8, align: 8, elements: !419, runtimeLang: DW_LANG_C_plus_plus)
!421 = !{ !423 }
!422 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "nested_exception", line: 1, size: 128, align: 64, elements: !421, runtimeLang: DW_LANG_C_plus_plus)
!423 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !422, name: "_M_ptr", line: 1, size: 64, align: 64, offset: 64, baseType: !15)
!424 = !{  }
!425 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "piecewise_construct_t", line: 1, size: 8, align: 8, elements: !424, runtimeLang: DW_LANG_C_plus_plus)
!426 = !{ !428 }
!427 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__nonesuch_no_braces", line: 1, size: 8, align: 8, elements: !426, runtimeLang: DW_LANG_C_plus_plus)
!428 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !427, name: "__nonesuch", line: 1, size: 8, align: 8, baseType: !420)
!429 = !{  }
!430 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "input_iterator_tag", line: 1, size: 8, align: 8, elements: !429, runtimeLang: DW_LANG_C_plus_plus)
!431 = !{  }
!432 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "output_iterator_tag", line: 1, size: 8, align: 8, elements: !431, runtimeLang: DW_LANG_C_plus_plus)
!433 = !{ !435 }
!434 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "forward_iterator_tag", line: 1, size: 8, align: 8, elements: !433, runtimeLang: DW_LANG_C_plus_plus)
!435 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !434, name: "input_iterator_tag", line: 1, size: 8, align: 8, baseType: !430)
!436 = !{ !438 }
!437 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "bidirectional_iterator_tag", line: 1, size: 8, align: 8, elements: !436, runtimeLang: DW_LANG_C_plus_plus)
!438 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !437, name: "forward_iterator_tag", line: 1, size: 8, align: 8, baseType: !434)
!439 = !{ !441 }
!440 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "random_access_iterator_tag", line: 1, size: 8, align: 8, elements: !439, runtimeLang: DW_LANG_C_plus_plus)
!441 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !440, name: "bidirectional_iterator_tag", line: 1, size: 8, align: 8, baseType: !437)
!442 = !{  }
!443 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "__undefined", line: 1, align: 8, elements: !442, runtimeLang: DW_LANG_C_plus_plus)
!444 = !{  }
!445 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "char_traits", line: 1, size: 8, align: 8, elements: !444, runtimeLang: DW_LANG_C_plus_plus)
!446 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !445, name: "char_type", line: 1, size: 8, align: 8, baseType: !43)
!447 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !445, name: "int_type", line: 1, size: 32, align: 32, baseType: !61)
!448 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !445, name: "pos_type", line: 1, align: 8, baseType: !112)
!449 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !445, name: "off_type", line: 1, size: 64, align: 64, baseType: !32)
!450 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !11, name: "__c_locale", line: 1, size: 64, align: 64, baseType: !253)
!451 = !{  }
!452 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__is_transparent", line: 1, align: 8, elements: !451, runtimeLang: DW_LANG_C_plus_plus)
!453 = !{  }
!454 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__allocator_traits_base", line: 1, size: 8, align: 8, elements: !453, runtimeLang: DW_LANG_C_plus_plus)
!455 = !{  }
!456 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "_Hash_impl", line: 1, size: 8, align: 8, elements: !455, runtimeLang: DW_LANG_C_plus_plus)
!457 = !{  }
!458 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "_Fnv_hash_impl", line: 1, size: 8, align: 8, elements: !457, runtimeLang: DW_LANG_C_plus_plus)
!459 = !{ !461, !462, !463, !464, !465, !466, !467, !468, !469, !470, !471, !472, !473, !474, !475, !476, !477, !478, !479, !480, !481, !482, !483, !484, !485, !486, !487, !488, !489, !490, !491, !492, !493, !494, !495, !496, !497, !498, !499, !500, !501, !502, !503, !504, !505, !506, !507, !508, !509, !510, !511, !512, !513, !514, !515, !516, !517, !518, !519, !520, !521, !522, !523, !524, !525, !526, !527, !528, !529, !530, !531, !532, !533, !534, !535, !536, !537, !538 }
!460 = !DICompositeType(tag: DW_TAG_enumeration_type, file: !3, scope: !11, name: "errc", line: 1, size: 32, align: 32, elements: !459, runtimeLang: DW_LANG_C_plus_plus)
!461 = !DIEnumerator(name: "_ZNSt4errc19wrong_protocol_typeE", value: 91)
!462 = !DIEnumerator(name: "_ZNSt4errc15value_too_largeE", value: 75)
!463 = !DIEnumerator(name: "_ZNSt4errc29too_many_symbolic_link_levelsE", value: 40)
!464 = !DIEnumerator(name: "_ZNSt4errc14too_many_linksE", value: 31)
!465 = !DIEnumerator(name: "_ZNSt4errc19too_many_files_openE", value: 24)
!466 = !DIEnumerator(name: "_ZNSt4errc29too_many_files_open_in_systemE", value: 23)
!467 = !DIEnumerator(name: "_ZNSt4errc9timed_outE", value: 110)
!468 = !DIEnumerator(name: "_ZNSt4errc14text_file_busyE", value: 26)
!469 = !DIEnumerator(name: "_ZNSt4errc14stream_timeoutE", value: 62)
!470 = !DIEnumerator(name: "_ZNSt4errc21state_not_recoverableE", value: 131)
!471 = !DIEnumerator(name: "_ZNSt4errc19result_out_of_rangeE", value: 34)
!472 = !DIEnumerator(name: "_ZNSt4errc30resource_unavailable_try_againE", value: 11)
!473 = !DIEnumerator(name: "_ZNSt4errc29resource_deadlock_would_occurE", value: 35)
!474 = !DIEnumerator(name: "_ZNSt4errc21read_only_file_systemE", value: 30)
!475 = !DIEnumerator(name: "_ZNSt4errc22protocol_not_supportedE", value: 93)
!476 = !DIEnumerator(name: "_ZNSt4errc14protocol_errorE", value: 71)
!477 = !DIEnumerator(name: "_ZNSt4errc17permission_deniedE", value: 13)
!478 = !DIEnumerator(name: "_ZNSt4errc10owner_deadE", value: 130)
!479 = !DIEnumerator(name: "_ZNSt4errc21operation_would_blockE", value: 11)
!480 = !DIEnumerator(name: "_ZNSt4errc23operation_not_supportedE", value: 95)
!481 = !DIEnumerator(name: "_ZNSt4errc23operation_not_permittedE", value: 1)
!482 = !DIEnumerator(name: "_ZNSt4errc21operation_in_progressE", value: 115)
!483 = !DIEnumerator(name: "_ZNSt4errc18operation_canceledE", value: 125)
!484 = !DIEnumerator(name: "_ZNSt4errc13not_supportedE", value: 95)
!485 = !DIEnumerator(name: "_ZNSt4errc17not_enough_memoryE", value: 12)
!486 = !DIEnumerator(name: "_ZNSt4errc13not_connectedE", value: 107)
!487 = !DIEnumerator(name: "_ZNSt4errc12not_a_streamE", value: 60)
!488 = !DIEnumerator(name: "_ZNSt4errc12not_a_socketE", value: 88)
!489 = !DIEnumerator(name: "_ZNSt4errc15not_a_directoryE", value: 20)
!490 = !DIEnumerator(name: "_ZNSt4errc15no_such_processE", value: 3)
!491 = !DIEnumerator(name: "_ZNSt4errc25no_such_file_or_directoryE", value: 2)
!492 = !DIEnumerator(name: "_ZNSt4errc14no_such_deviceE", value: 19)
!493 = !DIEnumerator(name: "_ZNSt4errc25no_such_device_or_addressE", value: 6)
!494 = !DIEnumerator(name: "_ZNSt4errc19no_stream_resourcesE", value: 63)
!495 = !DIEnumerator(name: "_ZNSt4errc18no_space_on_deviceE", value: 28)
!496 = !DIEnumerator(name: "_ZNSt4errc18no_protocol_optionE", value: 92)
!497 = !DIEnumerator(name: "_ZNSt4errc10no_messageE", value: 42)
!498 = !DIEnumerator(name: "_ZNSt4errc20no_message_availableE", value: 61)
!499 = !DIEnumerator(name: "_ZNSt4errc17no_lock_availableE", value: 37)
!500 = !DIEnumerator(name: "_ZNSt4errc7no_linkE", value: 67)
!501 = !DIEnumerator(name: "_ZNSt4errc16no_child_processE", value: 10)
!502 = !DIEnumerator(name: "_ZNSt4errc15no_buffer_spaceE", value: 105)
!503 = !DIEnumerator(name: "_ZNSt4errc19network_unreachableE", value: 101)
!504 = !DIEnumerator(name: "_ZNSt4errc13network_resetE", value: 102)
!505 = !DIEnumerator(name: "_ZNSt4errc12network_downE", value: 100)
!506 = !DIEnumerator(name: "_ZNSt4errc12message_sizeE", value: 90)
!507 = !DIEnumerator(name: "_ZNSt4errc14is_a_directoryE", value: 21)
!508 = !DIEnumerator(name: "_ZNSt4errc8io_errorE", value: 5)
!509 = !DIEnumerator(name: "_ZNSt4errc12invalid_seekE", value: 29)
!510 = !DIEnumerator(name: "_ZNSt4errc16invalid_argumentE", value: 22)
!511 = !DIEnumerator(name: "_ZNSt4errc11interruptedE", value: 4)
!512 = !DIEnumerator(name: "_ZNSt4errc34inappropriate_io_control_operationE", value: 25)
!513 = !DIEnumerator(name: "_ZNSt4errc21illegal_byte_sequenceE", value: 84)
!514 = !DIEnumerator(name: "_ZNSt4errc18identifier_removedE", value: 43)
!515 = !DIEnumerator(name: "_ZNSt4errc16host_unreachableE", value: 113)
!516 = !DIEnumerator(name: "_ZNSt4errc22function_not_supportedE", value: 38)
!517 = !DIEnumerator(name: "_ZNSt4errc17filename_too_longE", value: 36)
!518 = !DIEnumerator(name: "_ZNSt4errc14file_too_largeE", value: 27)
!519 = !DIEnumerator(name: "_ZNSt4errc11file_existsE", value: 17)
!520 = !DIEnumerator(name: "_ZNSt4errc23executable_format_errorE", value: 8)
!521 = !DIEnumerator(name: "_ZNSt4errc19directory_not_emptyE", value: 39)
!522 = !DIEnumerator(name: "_ZNSt4errc23device_or_resource_busyE", value: 16)
!523 = !DIEnumerator(name: "_ZNSt4errc28destination_address_requiredE", value: 89)
!524 = !DIEnumerator(name: "_ZNSt4errc17cross_device_linkE", value: 18)
!525 = !DIEnumerator(name: "_ZNSt4errc16connection_resetE", value: 104)
!526 = !DIEnumerator(name: "_ZNSt4errc18connection_refusedE", value: 111)
!527 = !DIEnumerator(name: "_ZNSt4errc30connection_already_in_progressE", value: 114)
!528 = !DIEnumerator(name: "_ZNSt4errc18connection_abortedE", value: 103)
!529 = !DIEnumerator(name: "_ZNSt4errc11broken_pipeE", value: 32)
!530 = !DIEnumerator(name: "_ZNSt4errc11bad_messageE", value: 74)
!531 = !DIEnumerator(name: "_ZNSt4errc19bad_file_descriptorE", value: 9)
!532 = !DIEnumerator(name: "_ZNSt4errc11bad_addressE", value: 14)
!533 = !DIEnumerator(name: "_ZNSt4errc22argument_out_of_domainE", value: 33)
!534 = !DIEnumerator(name: "_ZNSt4errc22argument_list_too_longE", value: 7)
!535 = !DIEnumerator(name: "_ZNSt4errc17already_connectedE", value: 106)
!536 = !DIEnumerator(name: "_ZNSt4errc21address_not_availableE", value: 99)
!537 = !DIEnumerator(name: "_ZNSt4errc14address_in_useE", value: 98)
!538 = !DIEnumerator(name: "_ZNSt4errc28address_family_not_supportedE", value: 97)
!539 = !{ !546 }
!540 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "__cow_string", line: 1, size: 64, align: 64, elements: !539, runtimeLang: DW_LANG_C_plus_plus)
!541 = !{ !543, !545 }
!542 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !540, line: 1, size: 64, align: 64, elements: !541, runtimeLang: DW_LANG_C_plus_plus)
!543 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !542, name: "_M_p", line: 1, size: 64, align: 64, baseType: !44)
!544 = !DICompositeType(tag: DW_TAG_array_type, size: 64, align: 8, baseType: !43, elements: !87)
!545 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !542, name: "_M_bytes", line: 1, size: 64, align: 8, baseType: !544)
!546 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !540, line: 1, size: 64, align: 64, baseType: !542)
!547 = !{ !549, !550 }
!548 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "logic_error", line: 1, size: 128, align: 64, elements: !547, runtimeLang: DW_LANG_C_plus_plus)
!549 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !548, name: "exception", line: 1, size: 64, align: 64, baseType: !347)
!550 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !548, name: "_M_msg", line: 1, size: 64, align: 64, offset: 64, baseType: !540)
!551 = !{ !553 }
!552 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "domain_error", line: 1, size: 128, align: 64, elements: !551, runtimeLang: DW_LANG_C_plus_plus)
!553 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !552, name: "logic_error", line: 1, size: 128, align: 64, baseType: !548)
!554 = !{ !556 }
!555 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "invalid_argument", line: 1, size: 128, align: 64, elements: !554, runtimeLang: DW_LANG_C_plus_plus)
!556 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !555, name: "logic_error", line: 1, size: 128, align: 64, baseType: !548)
!557 = !{ !559 }
!558 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "length_error", line: 1, size: 128, align: 64, elements: !557, runtimeLang: DW_LANG_C_plus_plus)
!559 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !558, name: "logic_error", line: 1, size: 128, align: 64, baseType: !548)
!560 = !{ !562 }
!561 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "out_of_range", line: 1, size: 128, align: 64, elements: !560, runtimeLang: DW_LANG_C_plus_plus)
!562 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !561, name: "logic_error", line: 1, size: 128, align: 64, baseType: !548)
!563 = !{ !565, !566 }
!564 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "runtime_error", line: 1, size: 128, align: 64, elements: !563, runtimeLang: DW_LANG_C_plus_plus)
!565 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !564, name: "exception", line: 1, size: 64, align: 64, baseType: !347)
!566 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !564, name: "_M_msg", line: 1, size: 64, align: 64, offset: 64, baseType: !540)
!567 = !{ !569 }
!568 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "range_error", line: 1, size: 128, align: 64, elements: !567, runtimeLang: DW_LANG_C_plus_plus)
!569 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !568, name: "runtime_error", line: 1, size: 128, align: 64, baseType: !564)
!570 = !{ !572 }
!571 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "overflow_error", line: 1, size: 128, align: 64, elements: !570, runtimeLang: DW_LANG_C_plus_plus)
!572 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !571, name: "runtime_error", line: 1, size: 128, align: 64, baseType: !564)
!573 = !{ !575 }
!574 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "underflow_error", line: 1, size: 128, align: 64, elements: !573, runtimeLang: DW_LANG_C_plus_plus)
!575 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !574, name: "runtime_error", line: 1, size: 128, align: 64, baseType: !564)
!576 = !{ !578, !583 }
!577 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "error_code", line: 1, size: 128, align: 64, elements: !576, runtimeLang: DW_LANG_C_plus_plus)
!578 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !577, name: "_M_value", line: 1, size: 32, align: 32, baseType: !61)
!579 = !{ !581 }
!580 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt3_V214error_categoryE", size: 64, align: 64, elements: !579)
!581 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !580, name: "__vptr", size: 64, align: 64, baseType: !126)
!582 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !580)
!583 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !577, name: "_M_cat", line: 1, size: 64, align: 64, offset: 64, baseType: !582)
!584 = !{ !586, !587 }
!585 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "error_condition", line: 1, size: 128, align: 64, elements: !584, runtimeLang: DW_LANG_C_plus_plus)
!586 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !585, name: "_M_value", line: 1, size: 32, align: 32, baseType: !61)
!587 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !585, name: "_M_cat", line: 1, size: 64, align: 64, offset: 64, baseType: !582)
!588 = !{ !590, !591 }
!589 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "system_error", line: 1, size: 256, align: 64, elements: !588, runtimeLang: DW_LANG_C_plus_plus)
!590 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !589, name: "runtime_error", line: 1, size: 128, align: 64, baseType: !564)
!591 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !589, name: "_M_code", line: 1, size: 128, align: 64, offset: 128, baseType: !577)
!592 = !{ !594, !595, !596, !597, !598, !599, !600, !601, !602 }
!593 = !DICompositeType(tag: DW_TAG_enumeration_type, file: !3, scope: !11, name: "_Ios_Openmode", line: 1, size: 32, align: 32, elements: !592, runtimeLang: DW_LANG_C_plus_plus)
!594 = !DIEnumerator(name: "_S_app", value: 1)
!595 = !DIEnumerator(name: "_S_ate", value: 2)
!596 = !DIEnumerator(name: "_S_bin", value: 4)
!597 = !DIEnumerator(name: "_S_in", value: 8)
!598 = !DIEnumerator(name: "_S_out", value: 16)
!599 = !DIEnumerator(name: "_S_trunc", value: 32)
!600 = !DIEnumerator(name: "_S_ios_openmode_end", value: 65536)
!601 = !DIEnumerator(name: "_S_ios_openmode_max", value: 2147483647)
!602 = !DIEnumerator(name: "_S_ios_openmode_min", value: -2147483648)
!603 = !{ !605, !606, !607, !608 }
!604 = !DICompositeType(tag: DW_TAG_enumeration_type, file: !3, scope: !11, name: "_Ios_Seekdir", line: 1, size: 32, align: 32, elements: !603, runtimeLang: DW_LANG_C_plus_plus)
!605 = !DIEnumerator(name: "_S_beg", value: 0)
!606 = !DIEnumerator(name: "_S_cur", value: 1)
!607 = !DIEnumerator(name: "_S_end", value: 2)
!608 = !DIEnumerator(name: "_S_ios_seekdir_end", value: 65536)
!609 = !{ !611 }
!610 = !DICompositeType(tag: DW_TAG_enumeration_type, file: !3, scope: !11, name: "io_errc", line: 1, size: 32, align: 32, elements: !609, runtimeLang: DW_LANG_C_plus_plus)
!611 = !DIEnumerator(name: "_ZNSt7io_errc6streamE", value: 1)
!612 = !{  }
!613 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !11, name: "ctype_base", line: 1, size: 8, align: 8, elements: !612, runtimeLang: DW_LANG_C_plus_plus)
!614 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !613, name: "__to_type", line: 1, size: 64, align: 64, baseType: !62)
!615 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !613, name: "mask", line: 1, size: 16, align: 16, baseType: !79)
!616 = !{  }
!617 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !11, name: "__num_base", line: 1, size: 8, align: 8, elements: !616, runtimeLang: DW_LANG_C_plus_plus)
!618 = !DINamespace(scope: !10, name: "__cxxabiv1")
!619 = !{  }
!620 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !618, name: "__cxa_refcounted_exception", line: 1, align: 8, elements: !619, runtimeLang: DW_LANG_C_plus_plus)
!621 = !{  }
!622 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !618, name: "__class_type_info", line: 1, align: 8, elements: !621, runtimeLang: DW_LANG_C_plus_plus)
!623 = !{  }
!624 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !618, name: "__forced_unwind", line: 1, size: 64, align: 64, elements: !623, runtimeLang: DW_LANG_C_plus_plus)
!625 = !DINamespace(scope: !10, name: "__gnu_cxx")
!626 = !{  }
!627 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, line: 1, size: 128, align: 64, elements: !626, runtimeLang: DW_LANG_C_plus_plus)
!628 = !{  }
!629 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__class_type_info", line: 1, size: 128, align: 64, elements: !628, runtimeLang: DW_LANG_C_plus_plus)
!630 = !{  }
!631 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__EDG_type_info", line: 1, size: 128, align: 64, elements: !630, runtimeLang: DW_LANG_C_plus_plus)
!632 = !{  }
!633 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pbase_type_info", line: 1, size: 256, align: 64, elements: !632, runtimeLang: DW_LANG_C_plus_plus)
!634 = !{  }
!635 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pointer_to_member_type_info", line: 1, size: 320, align: 64, elements: !634, runtimeLang: DW_LANG_C_plus_plus)
!636 = !{  }
!637 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pointer_type_info", line: 1, size: 256, align: 64, elements: !636, runtimeLang: DW_LANG_C_plus_plus)
!638 = !{  }
!639 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__vmi_class_type_info", line: 1, size: 448, align: 64, elements: !638, runtimeLang: DW_LANG_C_plus_plus)
!640 = !{  }
!641 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__si_class_type_info", line: 1, size: 192, align: 64, elements: !640, runtimeLang: DW_LANG_C_plus_plus)
!642 = !{  }
!643 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__function_type_info", line: 1, size: 128, align: 64, elements: !642, runtimeLang: DW_LANG_C_plus_plus)
!644 = !{  }
!645 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__array_type_info", line: 1, size: 128, align: 64, elements: !644, runtimeLang: DW_LANG_C_plus_plus)
!646 = !{  }
!647 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__enum_type_info", line: 1, size: 128, align: 64, elements: !646, runtimeLang: DW_LANG_C_plus_plus)
!648 = !{  }
!649 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__fundamental_type_info", line: 1, size: 128, align: 64, elements: !648, runtimeLang: DW_LANG_C_plus_plus)
!650 = !DIBasicType(tag: DW_TAG_base_type, name: "__float128", size: 128, align: 128, encoding: DW_ATE_float)
!651 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__float128", line: 1, size: 128, align: 128, baseType: !650)
!652 = !{ !654, !655, !656, !657 }
!653 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__va_list_tag", line: 1, size: 192, align: 64, elements: !652, runtimeLang: DW_LANG_C_plus_plus)
!654 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !653, name: "gp_offset", line: 1, size: 32, align: 32, baseType: !97)
!655 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !653, name: "fp_offset", line: 1, size: 32, align: 32, offset: 32, baseType: !97)
!656 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !653, name: "overflow_arg_area", line: 1, size: 64, align: 64, offset: 64, baseType: !44)
!657 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !653, name: "reg_save_area", line: 1, size: 64, align: 64, offset: 128, baseType: !44)
!658 = !DICompositeType(tag: DW_TAG_array_type, size: 192, align: 64, baseType: !653, elements: !377)
!659 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__pgi_va_list", line: 1, size: 192, align: 64, baseType: !658)
!660 = !DIBasicType(tag: DW_TAG_base_type, name: "unsigned char", size: 8, align: 8, encoding: DW_ATE_unsigned_char)
!661 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__u_char", line: 1, size: 8, align: 8, baseType: !660)
!662 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__u_short", line: 1, size: 16, align: 16, baseType: !79)
!663 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__u_int", line: 1, size: 32, align: 32, baseType: !97)
!664 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__u_long", line: 1, size: 64, align: 64, baseType: !30)
!665 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__int8_t", line: 1, size: 8, align: 8, baseType: !43)
!666 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uint8_t", line: 1, size: 8, align: 8, baseType: !660)
!667 = !DIBasicType(tag: DW_TAG_base_type, name: "short", size: 16, align: 16, encoding: DW_ATE_signed)
!668 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__int16_t", line: 1, size: 16, align: 16, baseType: !667)
!669 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uint16_t", line: 1, size: 16, align: 16, baseType: !79)
!670 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__int32_t", line: 1, size: 32, align: 32, baseType: !61)
!671 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uint32_t", line: 1, size: 32, align: 32, baseType: !97)
!672 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__int64_t", line: 1, size: 64, align: 64, baseType: !32)
!673 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uint64_t", line: 1, size: 64, align: 64, baseType: !30)
!674 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__int_least8_t", line: 1, size: 8, align: 8, baseType: !43)
!675 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uint_least8_t", line: 1, size: 8, align: 8, baseType: !660)
!676 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__int_least16_t", line: 1, size: 16, align: 16, baseType: !667)
!677 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uint_least16_t", line: 1, size: 16, align: 16, baseType: !79)
!678 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__int_least32_t", line: 1, size: 32, align: 32, baseType: !61)
!679 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uint_least32_t", line: 1, size: 32, align: 32, baseType: !97)
!680 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__int_least64_t", line: 1, size: 64, align: 64, baseType: !32)
!681 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uint_least64_t", line: 1, size: 64, align: 64, baseType: !30)
!682 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__quad_t", line: 1, size: 64, align: 64, baseType: !32)
!683 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__u_quad_t", line: 1, size: 64, align: 64, baseType: !30)
!684 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__intmax_t", line: 1, size: 64, align: 64, baseType: !32)
!685 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uintmax_t", line: 1, size: 64, align: 64, baseType: !30)
!686 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__dev_t", line: 1, size: 64, align: 64, baseType: !30)
!687 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__uid_t", line: 1, size: 32, align: 32, baseType: !97)
!688 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gid_t", line: 1, size: 32, align: 32, baseType: !97)
!689 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__ino_t", line: 1, size: 64, align: 64, baseType: !30)
!690 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__ino64_t", line: 1, size: 64, align: 64, baseType: !30)
!691 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__mode_t", line: 1, size: 32, align: 32, baseType: !97)
!692 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__nlink_t", line: 1, size: 64, align: 64, baseType: !30)
!693 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__off_t", line: 1, size: 64, align: 64, baseType: !32)
!694 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__pid_t", line: 1, size: 32, align: 32, baseType: !61)
!695 = !{ !700 }
!696 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__fsid_t", line: 1, size: 64, align: 32, elements: !695, runtimeLang: DW_LANG_C_plus_plus)
!697 = !DISubrange(count: 2)
!698 = !{ !697 }
!699 = !DICompositeType(tag: DW_TAG_array_type, size: 64, align: 32, baseType: !61, elements: !698)
!700 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !696, name: "__val", line: 1, size: 64, align: 32, baseType: !699)
!701 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fsid_t", line: 1, size: 64, align: 32, baseType: !696)
!702 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__clock_t", line: 1, size: 64, align: 64, baseType: !32)
!703 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__rlim_t", line: 1, size: 64, align: 64, baseType: !30)
!704 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__rlim64_t", line: 1, size: 64, align: 64, baseType: !30)
!705 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__id_t", line: 1, size: 32, align: 32, baseType: !97)
!706 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__time_t", line: 1, size: 64, align: 64, baseType: !32)
!707 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__useconds_t", line: 1, size: 32, align: 32, baseType: !97)
!708 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__suseconds_t", line: 1, size: 64, align: 64, baseType: !32)
!709 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__daddr_t", line: 1, size: 32, align: 32, baseType: !61)
!710 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__key_t", line: 1, size: 32, align: 32, baseType: !61)
!711 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__clockid_t", line: 1, size: 32, align: 32, baseType: !61)
!712 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__timer_t", line: 1, size: 64, align: 64, baseType: !17)
!713 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__blksize_t", line: 1, size: 64, align: 64, baseType: !32)
!714 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__blkcnt_t", line: 1, size: 64, align: 64, baseType: !32)
!715 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__blkcnt64_t", line: 1, size: 64, align: 64, baseType: !32)
!716 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fsblkcnt_t", line: 1, size: 64, align: 64, baseType: !30)
!717 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fsblkcnt64_t", line: 1, size: 64, align: 64, baseType: !30)
!718 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fsfilcnt_t", line: 1, size: 64, align: 64, baseType: !30)
!719 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fsfilcnt64_t", line: 1, size: 64, align: 64, baseType: !30)
!720 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fsword_t", line: 1, size: 64, align: 64, baseType: !32)
!721 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__syscall_slong_t", line: 1, size: 64, align: 64, baseType: !32)
!722 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__syscall_ulong_t", line: 1, size: 64, align: 64, baseType: !30)
!723 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__loff_t", line: 1, size: 64, align: 64, baseType: !32)
!724 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__caddr_t", line: 1, size: 64, align: 64, baseType: !44)
!725 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__intptr_t", line: 1, size: 64, align: 64, baseType: !32)
!726 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__socklen_t", line: 1, size: 32, align: 32, baseType: !97)
!727 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__sig_atomic_t", line: 1, size: 32, align: 32, baseType: !61)
!728 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_Float128", line: 1, size: 128, align: 128, baseType: !650)
!729 = !DIBasicType(tag: DW_TAG_base_type, name: "float", size: 32, align: 32, encoding: DW_ATE_float)
!730 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_Float32", line: 1, size: 32, align: 32, baseType: !729)
!731 = !DIBasicType(tag: DW_TAG_base_type, name: "double", size: 64, align: 64, encoding: DW_ATE_float)
!732 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_Float64", line: 1, size: 64, align: 64, baseType: !731)
!733 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_Float32x", line: 1, size: 64, align: 64, baseType: !731)
!734 = !DIBasicType(tag: DW_TAG_base_type, name: "80-bit extended precision", size: 128, align: 128, encoding: DW_ATE_signed)
!735 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_Float64x", line: 1, size: 128, align: 128, baseType: !734)
!736 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "double_t", line: 1, size: 64, align: 64, baseType: !731)
!737 = !{ !739, !740, !741 }
!738 = !DICompositeType(tag: DW_TAG_enumeration_type, file: !3, scope: !10, name: "idtype_t", line: 1, size: 32, align: 32, elements: !737, runtimeLang: DW_LANG_C_plus_plus)
!739 = !DIEnumerator(name: "P_ALL", value: 0)
!740 = !DIEnumerator(name: "P_PID", value: 1)
!741 = !DIEnumerator(name: "P_PGID", value: 2)
!742 = !{ !744, !745 }
!743 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "div_t", line: 1, size: 64, align: 32, elements: !742, runtimeLang: DW_LANG_C_plus_plus)
!744 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !743, name: "quot", line: 1, size: 32, align: 32, baseType: !61)
!745 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !743, name: "rem", line: 1, size: 32, align: 32, offset: 32, baseType: !61)
!746 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "div_t", line: 1, size: 64, align: 32, baseType: !743)
!747 = !{ !749, !750 }
!748 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "ldiv_t", line: 1, size: 128, align: 64, elements: !747, runtimeLang: DW_LANG_C_plus_plus)
!749 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !748, name: "quot", line: 1, size: 64, align: 64, baseType: !32)
!750 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !748, name: "rem", line: 1, size: 64, align: 64, offset: 64, baseType: !32)
!751 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "ldiv_t", line: 1, size: 128, align: 64, baseType: !748)
!752 = !{ !755, !756 }
!753 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "lldiv_t", line: 1, size: 128, align: 64, elements: !752, runtimeLang: DW_LANG_C_plus_plus)
!754 = !DIBasicType(tag: DW_TAG_base_type, name: "long long", size: 64, align: 64, encoding: DW_ATE_signed)
!755 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !753, name: "quot", line: 1, size: 64, align: 64, baseType: !754)
!756 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !753, name: "rem", line: 1, size: 64, align: 64, offset: 64, baseType: !754)
!757 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "lldiv_t", line: 1, size: 128, align: 64, baseType: !753)
!758 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__locale_t", line: 1, size: 64, align: 64, baseType: !253)
!759 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "locale_t", line: 1, size: 64, align: 64, baseType: !253)
!760 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pid_t", line: 1, size: 32, align: 32, baseType: !61)
!761 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "clock_t", line: 1, size: 64, align: 64, baseType: !32)
!762 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "clockid_t", line: 1, size: 32, align: 32, baseType: !61)
!763 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "time_t", line: 1, size: 64, align: 64, baseType: !32)
!764 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "timer_t", line: 1, size: 64, align: 64, baseType: !17)
!765 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "ulong", line: 1, size: 64, align: 64, baseType: !30)
!766 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "uint", line: 1, size: 32, align: 32, baseType: !97)
!767 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "int32_t", line: 1, size: 32, align: 32, baseType: !61)
!768 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "register_t", line: 1, size: 64, align: 64, baseType: !32)
!769 = !{ !771 }
!770 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__sigset_t", line: 1, size: 1024, align: 64, elements: !769, runtimeLang: DW_LANG_C_plus_plus)
!771 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !770, name: "__val", line: 1, size: 1024, align: 64, baseType: !329)
!772 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__sigset_t", line: 1, size: 1024, align: 64, baseType: !770)
!773 = !{ !775, !776 }
!774 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "timeval", line: 1, size: 128, align: 64, elements: !773, runtimeLang: DW_LANG_C_plus_plus)
!775 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !774, name: "tv_sec", line: 1, size: 64, align: 64, baseType: !32)
!776 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !774, name: "tv_usec", line: 1, size: 64, align: 64, offset: 64, baseType: !32)
!777 = !{ !779, !780 }
!778 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "timespec", line: 1, size: 128, align: 64, elements: !777, runtimeLang: DW_LANG_C_plus_plus)
!779 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !778, name: "tv_sec", line: 1, size: 64, align: 64, baseType: !32)
!780 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !778, name: "tv_nsec", line: 1, size: 64, align: 64, offset: 64, baseType: !32)
!781 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fd_mask", line: 1, size: 64, align: 64, baseType: !32)
!782 = !{ !785 }
!783 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "fd_set", line: 1, size: 1024, align: 64, elements: !782, runtimeLang: DW_LANG_C_plus_plus)
!784 = !DICompositeType(tag: DW_TAG_array_type, size: 1024, align: 64, baseType: !32, elements: !51)
!785 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !783, name: "fds_bits", line: 1, size: 1024, align: 64, baseType: !784)
!786 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "fd_set", line: 1, size: 1024, align: 64, baseType: !783)
!787 = !{ !790, !791 }
!788 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pthread_internal_list", line: 1, size: 128, align: 64, elements: !787, runtimeLang: DW_LANG_C_plus_plus)
!789 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !788)
!790 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !788, name: "__prev", line: 1, size: 64, align: 64, baseType: !789)
!791 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !788, name: "__next", line: 1, size: 64, align: 64, offset: 64, baseType: !789)
!792 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__pthread_list_t", line: 1, size: 128, align: 64, baseType: !788)
!793 = !{ !796 }
!794 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pthread_internal_slist", line: 1, size: 64, align: 64, elements: !793, runtimeLang: DW_LANG_C_plus_plus)
!795 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !794)
!796 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !794, name: "__next", line: 1, size: 64, align: 64, baseType: !795)
!797 = !{ !799, !800, !801, !802, !803, !804, !805, !806 }
!798 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pthread_mutex_s", line: 1, size: 320, align: 64, elements: !797, runtimeLang: DW_LANG_C_plus_plus)
!799 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !798, name: "__lock", line: 1, size: 32, align: 32, baseType: !61)
!800 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !798, name: "__count", line: 1, size: 32, align: 32, offset: 32, baseType: !97)
!801 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !798, name: "__owner", line: 1, size: 32, align: 32, offset: 64, baseType: !61)
!802 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !798, name: "__nusers", line: 1, size: 32, align: 32, offset: 96, baseType: !97)
!803 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !798, name: "__kind", line: 1, size: 32, align: 32, offset: 128, baseType: !61)
!804 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !798, name: "__spins", line: 1, size: 16, align: 16, offset: 160, baseType: !667)
!805 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !798, name: "__elision", line: 1, size: 16, align: 16, offset: 176, baseType: !667)
!806 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !798, name: "__list", line: 1, size: 128, align: 64, offset: 192, baseType: !788)
!807 = !{ !809, !810, !811, !812, !813, !814, !815, !816, !817, !821, !822, !823 }
!808 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pthread_rwlock_arch_t", line: 1, size: 448, align: 64, elements: !807, runtimeLang: DW_LANG_C_plus_plus)
!809 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__readers", line: 1, size: 32, align: 32, baseType: !97)
!810 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__writers", line: 1, size: 32, align: 32, offset: 32, baseType: !97)
!811 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__wrphase_futex", line: 1, size: 32, align: 32, offset: 64, baseType: !97)
!812 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__writers_futex", line: 1, size: 32, align: 32, offset: 96, baseType: !97)
!813 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__pad3", line: 1, size: 32, align: 32, offset: 128, baseType: !97)
!814 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__pad4", line: 1, size: 32, align: 32, offset: 160, baseType: !97)
!815 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__cur_writer", line: 1, size: 32, align: 32, offset: 192, baseType: !61)
!816 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__shared", line: 1, size: 32, align: 32, offset: 224, baseType: !61)
!817 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__rwelision", line: 1, size: 8, align: 8, offset: 256, baseType: !43)
!818 = !DISubrange(count: 7)
!819 = !{ !818 }
!820 = !DICompositeType(tag: DW_TAG_array_type, size: 56, align: 8, baseType: !660, elements: !819)
!821 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__pad1", line: 1, size: 56, align: 8, offset: 264, baseType: !820)
!822 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__pad2", line: 1, size: 64, align: 64, offset: 320, baseType: !30)
!823 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !808, name: "__flags", line: 1, size: 32, align: 32, offset: 384, baseType: !97)
!824 = !{ !843, !843, !845, !846, !847, !848, !849 }
!825 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pthread_cond_s", line: 1, size: 384, align: 64, elements: !824, runtimeLang: DW_LANG_C_plus_plus)
!826 = !{ !833, !834 }
!827 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !825, line: 1, size: 64, align: 64, elements: !826, runtimeLang: DW_LANG_C_plus_plus)
!828 = !{ !830, !831 }
!829 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !827, name: "_ZN16__pthread_cond_s4__C2Ut_E", line: 1, size: 64, align: 32, elements: !828, runtimeLang: DW_LANG_C_plus_plus)
!830 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !829, name: "__low", line: 1, size: 32, align: 32, baseType: !97)
!831 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !829, name: "__high", line: 1, size: 32, align: 32, offset: 32, baseType: !97)
!832 = !DIBasicType(tag: DW_TAG_base_type, name: "unsigned long long", size: 64, align: 64, encoding: DW_ATE_unsigned)
!833 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !827, name: "__wseq", line: 1, size: 64, align: 64, baseType: !832)
!834 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !827, name: "__wseq32", line: 1, size: 64, align: 32, baseType: !829)
!835 = !{ !841, !842 }
!836 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !825, line: 1, size: 64, align: 64, elements: !835, runtimeLang: DW_LANG_C_plus_plus)
!837 = !{ !839, !840 }
!838 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !836, name: "_ZN16__pthread_cond_s4__C3Ut_E", line: 1, size: 64, align: 32, elements: !837, runtimeLang: DW_LANG_C_plus_plus)
!839 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !838, name: "__low", line: 1, size: 32, align: 32, baseType: !97)
!840 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !838, name: "__high", line: 1, size: 32, align: 32, offset: 32, baseType: !97)
!841 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !836, name: "__g1_start", line: 1, size: 64, align: 64, baseType: !832)
!842 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !836, name: "__g1_start32", line: 1, size: 64, align: 32, baseType: !838)
!843 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !825, line: 1, size: 64, align: 64, baseType: !827)
!844 = !DICompositeType(tag: DW_TAG_array_type, size: 64, align: 32, baseType: !97, elements: !698)
!845 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !825, name: "__g_refs", line: 1, size: 64, align: 32, offset: 128, baseType: !844)
!846 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !825, name: "__g_size", line: 1, size: 64, align: 32, offset: 192, baseType: !844)
!847 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !825, name: "__g1_orig_size", line: 1, size: 32, align: 32, offset: 256, baseType: !97)
!848 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !825, name: "__wrefs", line: 1, size: 32, align: 32, offset: 288, baseType: !97)
!849 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !825, name: "__g_signals", line: 1, size: 64, align: 32, offset: 320, baseType: !844)
!850 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_t", line: 1, size: 64, align: 64, baseType: !30)
!851 = !{ !854, !855 }
!852 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_mutexattr_t", line: 1, size: 32, align: 32, elements: !851, runtimeLang: DW_LANG_C_plus_plus)
!853 = !DICompositeType(tag: DW_TAG_array_type, size: 32, align: 8, baseType: !43, elements: !69)
!854 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !852, name: "__size", line: 1, size: 32, align: 8, baseType: !853)
!855 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !852, name: "__align", line: 1, size: 32, align: 32, baseType: !61)
!856 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_mutexattr_t", line: 1, size: 32, align: 32, baseType: !852)
!857 = !{ !859, !860 }
!858 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_condattr_t", line: 1, size: 32, align: 32, elements: !857, runtimeLang: DW_LANG_C_plus_plus)
!859 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !858, name: "__size", line: 1, size: 32, align: 8, baseType: !853)
!860 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !858, name: "__align", line: 1, size: 32, align: 32, baseType: !61)
!861 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_condattr_t", line: 1, size: 32, align: 32, baseType: !858)
!862 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_key_t", line: 1, size: 32, align: 32, baseType: !97)
!863 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_once_t", line: 1, size: 32, align: 32, baseType: !61)
!864 = !{ !869, !870 }
!865 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_attr_t", line: 1, size: 448, align: 64, elements: !864, runtimeLang: DW_LANG_C_plus_plus)
!866 = !DISubrange(count: 56)
!867 = !{ !866 }
!868 = !DICompositeType(tag: DW_TAG_array_type, size: 448, align: 8, baseType: !43, elements: !867)
!869 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !865, name: "__size", line: 1, size: 448, align: 8, baseType: !868)
!870 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !865, name: "__align", line: 1, size: 64, align: 64, baseType: !32)
!871 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_attr_t", line: 1, size: 448, align: 64, baseType: !865)
!872 = !{ !874, !878, !879 }
!873 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_mutex_t", line: 1, size: 320, align: 64, elements: !872, runtimeLang: DW_LANG_C_plus_plus)
!874 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !873, name: "__data", line: 1, size: 320, align: 64, baseType: !798)
!875 = !DISubrange(count: 40)
!876 = !{ !875 }
!877 = !DICompositeType(tag: DW_TAG_array_type, size: 320, align: 8, baseType: !43, elements: !876)
!878 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !873, name: "__size", line: 1, size: 320, align: 8, baseType: !877)
!879 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !873, name: "__align", line: 1, size: 64, align: 64, baseType: !32)
!880 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_mutex_t", line: 1, size: 320, align: 64, baseType: !873)
!881 = !{ !883, !887, !888 }
!882 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_cond_t", line: 1, size: 384, align: 64, elements: !881, runtimeLang: DW_LANG_C_plus_plus)
!883 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !882, name: "__data", line: 1, size: 384, align: 64, baseType: !825)
!884 = !DISubrange(count: 48)
!885 = !{ !884 }
!886 = !DICompositeType(tag: DW_TAG_array_type, size: 384, align: 8, baseType: !43, elements: !885)
!887 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !882, name: "__size", line: 1, size: 384, align: 8, baseType: !886)
!888 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !882, name: "__align", line: 1, size: 64, align: 64, baseType: !754)
!889 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_cond_t", line: 1, size: 384, align: 64, baseType: !882)
!890 = !{ !892, !893, !894 }
!891 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_rwlock_t", line: 1, size: 448, align: 64, elements: !890, runtimeLang: DW_LANG_C_plus_plus)
!892 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !891, name: "__data", line: 1, size: 448, align: 64, baseType: !808)
!893 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !891, name: "__size", line: 1, size: 448, align: 8, baseType: !868)
!894 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !891, name: "__align", line: 1, size: 64, align: 64, baseType: !32)
!895 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_rwlock_t", line: 1, size: 448, align: 64, baseType: !891)
!896 = !{ !898, !899 }
!897 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_rwlockattr_t", line: 1, size: 64, align: 64, elements: !896, runtimeLang: DW_LANG_C_plus_plus)
!898 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !897, name: "__size", line: 1, size: 64, align: 8, baseType: !544)
!899 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !897, name: "__align", line: 1, size: 64, align: 64, baseType: !32)
!900 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_rwlockattr_t", line: 1, size: 64, align: 64, baseType: !897)
!901 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_spinlock_t", line: 1, size: 32, align: 32, baseType: !61)
!902 = !{ !907, !908 }
!903 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_barrier_t", line: 1, size: 256, align: 64, elements: !902, runtimeLang: DW_LANG_C_plus_plus)
!904 = !DISubrange(count: 32)
!905 = !{ !904 }
!906 = !DICompositeType(tag: DW_TAG_array_type, size: 256, align: 8, baseType: !43, elements: !905)
!907 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !903, name: "__size", line: 1, size: 256, align: 8, baseType: !906)
!908 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !903, name: "__align", line: 1, size: 64, align: 64, baseType: !32)
!909 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_barrier_t", line: 1, size: 256, align: 64, baseType: !903)
!910 = !{ !912, !913 }
!911 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !10, name: "pthread_barrierattr_t", line: 1, size: 32, align: 32, elements: !910, runtimeLang: DW_LANG_C_plus_plus)
!912 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !911, name: "__size", line: 1, size: 32, align: 8, baseType: !853)
!913 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !911, name: "__align", line: 1, size: 32, align: 32, baseType: !61)
!914 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "pthread_barrierattr_t", line: 1, size: 32, align: 32, baseType: !911)
!915 = !{ !917, !918, !919, !920, !921, !922, !923 }
!916 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "random_data", line: 1, size: 384, align: 64, elements: !915, runtimeLang: DW_LANG_C_plus_plus)
!917 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !916, name: "fptr", line: 1, size: 64, align: 64, baseType: !62)
!918 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !916, name: "rptr", line: 1, size: 64, align: 64, offset: 64, baseType: !62)
!919 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !916, name: "state", line: 1, size: 64, align: 64, offset: 128, baseType: !62)
!920 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !916, name: "rand_type", line: 1, size: 32, align: 32, offset: 192, baseType: !61)
!921 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !916, name: "rand_deg", line: 1, size: 32, align: 32, offset: 224, baseType: !61)
!922 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !916, name: "rand_sep", line: 1, size: 32, align: 32, offset: 256, baseType: !61)
!923 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !916, name: "end_ptr", line: 1, size: 64, align: 64, offset: 320, baseType: !62)
!924 = !{ !929, !930, !931, !932, !933 }
!925 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "drand48_data", line: 1, size: 192, align: 64, elements: !924, runtimeLang: DW_LANG_C_plus_plus)
!926 = !DISubrange(count: 3)
!927 = !{ !926 }
!928 = !DICompositeType(tag: DW_TAG_array_type, size: 48, align: 16, baseType: !79, elements: !927)
!929 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !925, name: "__x", line: 1, size: 48, align: 16, baseType: !928)
!930 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !925, name: "__old_x", line: 1, size: 48, align: 16, offset: 48, baseType: !928)
!931 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !925, name: "__c", line: 1, size: 16, align: 16, offset: 96, baseType: !79)
!932 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !925, name: "__init", line: 1, size: 16, align: 16, offset: 112, baseType: !79)
!933 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !925, name: "__a", line: 1, size: 64, align: 64, offset: 128, baseType: !832)
!934 = !{ !61, !17, !17 }
!935 = !DISubroutineType(types: !934)
!936 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !935)
!937 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__compar_fn_t", line: 1, size: 64, align: 64, baseType: !936)
!938 = !{ !61, !17, !17, !17 }
!939 = !DISubroutineType(types: !938)
!940 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !939)
!941 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__compar_d_fn_t", line: 1, size: 64, align: 64, baseType: !940)
!942 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gnuc_va_list", line: 1, size: 192, align: 64, baseType: !658)
!943 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "wint_t", line: 1, size: 32, align: 32, baseType: !97)
!944 = !{ !950, !951 }
!945 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__mbstate_t", line: 1, size: 64, align: 32, elements: !944, runtimeLang: DW_LANG_C_plus_plus)
!946 = !{ !948, !949 }
!947 = !DICompositeType(tag: DW_TAG_union_type, file: !3, scope: !945, name: "_ZN11__mbstate_tUt_E", line: 1, size: 32, align: 32, elements: !946, runtimeLang: DW_LANG_C_plus_plus)
!948 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !947, name: "__wch", line: 1, size: 32, align: 32, baseType: !97)
!949 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !947, name: "__wchb", line: 1, size: 32, align: 8, baseType: !853)
!950 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !945, name: "__count", line: 1, size: 32, align: 32, baseType: !61)
!951 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !945, name: "__value", line: 1, size: 32, align: 32, offset: 32, baseType: !947)
!952 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__mbstate_t", line: 1, size: 64, align: 32, baseType: !945)
!953 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "mbstate_t", line: 1, size: 64, align: 32, baseType: !945)
!954 = !{ !956, !957, !958, !959, !960, !961, !962, !963, !964, !965, !966, !967, !971, !973, !974, !975, !976, !977, !978, !979, !980, !981, !985, !989, !990, !991, !992, !993, !997 }
!955 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_IO_FILE", line: 1, size: 1728, align: 64, elements: !954, runtimeLang: DW_LANG_C_plus_plus)
!956 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_flags", line: 1, size: 32, align: 32, baseType: !61)
!957 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_read_ptr", line: 1, size: 64, align: 64, offset: 64, baseType: !44)
!958 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_read_end", line: 1, size: 64, align: 64, offset: 128, baseType: !44)
!959 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_read_base", line: 1, size: 64, align: 64, offset: 192, baseType: !44)
!960 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_write_base", line: 1, size: 64, align: 64, offset: 256, baseType: !44)
!961 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_write_ptr", line: 1, size: 64, align: 64, offset: 320, baseType: !44)
!962 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_write_end", line: 1, size: 64, align: 64, offset: 384, baseType: !44)
!963 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_buf_base", line: 1, size: 64, align: 64, offset: 448, baseType: !44)
!964 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_buf_end", line: 1, size: 64, align: 64, offset: 512, baseType: !44)
!965 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_save_base", line: 1, size: 64, align: 64, offset: 576, baseType: !44)
!966 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_backup_base", line: 1, size: 64, align: 64, offset: 640, baseType: !44)
!967 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_IO_save_end", line: 1, size: 64, align: 64, offset: 704, baseType: !44)
!968 = !{  }
!969 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_IO_marker", align: 8, elements: !968)
!970 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !969)
!971 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_markers", line: 1, size: 64, align: 64, offset: 768, baseType: !970)
!972 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !955)
!973 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_chain", line: 1, size: 64, align: 64, offset: 832, baseType: !972)
!974 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_fileno", line: 1, size: 32, align: 32, offset: 896, baseType: !61)
!975 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_flags2", line: 1, size: 32, align: 32, offset: 928, baseType: !61)
!976 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_old_offset", line: 1, size: 64, align: 64, offset: 960, baseType: !32)
!977 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_cur_column", line: 1, size: 16, align: 16, offset: 1024, baseType: !79)
!978 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_vtable_offset", line: 1, size: 8, align: 8, offset: 1040, baseType: !43)
!979 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_shortbuf", line: 1, size: 8, align: 8, offset: 1048, baseType: !378)
!980 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_lock", line: 1, size: 64, align: 64, offset: 1088, baseType: !17)
!981 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_offset", line: 1, size: 64, align: 64, offset: 1152, baseType: !32)
!982 = !{  }
!983 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_IO_codecvt", align: 8, elements: !982)
!984 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !983)
!985 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_codecvt", line: 1, size: 64, align: 64, offset: 1216, baseType: !984)
!986 = !{  }
!987 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_IO_wide_data", align: 8, elements: !986)
!988 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !987)
!989 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_wide_data", line: 1, size: 64, align: 64, offset: 1280, baseType: !988)
!990 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_freeres_list", line: 1, size: 64, align: 64, offset: 1344, baseType: !972)
!991 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_freeres_buf", line: 1, size: 64, align: 64, offset: 1408, baseType: !17)
!992 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "__pad5", line: 1, size: 64, align: 64, offset: 1472, baseType: !30)
!993 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_mode", line: 1, size: 32, align: 32, offset: 1536, baseType: !61)
!994 = !DISubrange(count: 20)
!995 = !{ !994 }
!996 = !DICompositeType(tag: DW_TAG_array_type, size: 160, align: 8, baseType: !43, elements: !995)
!997 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !955, name: "_unused2", line: 1, size: 160, align: 8, offset: 1568, baseType: !996)
!998 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__FILE", line: 1, size: 1728, align: 64, baseType: !955)
!999 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "FILE", line: 1, size: 1728, align: 64, baseType: !955)
!1000 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "ptrdiff_t", line: 1, size: 64, align: 64, baseType: !32)
!1001 = !{ !1003, !1004 }
!1002 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "max_align_t", line: 1, size: 256, align: 128, elements: !1001, runtimeLang: DW_LANG_C_plus_plus)
!1003 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1002, name: "__max_align_ll", line: 1, size: 64, align: 64, baseType: !754)
!1004 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1002, name: "__max_align_ld", line: 1, size: 128, align: 128, offset: 128, baseType: !734)
!1005 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "uint_least16_t", line: 1, size: 16, align: 16, baseType: !79)
!1006 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "uint_least32_t", line: 1, size: 32, align: 32, baseType: !97)
!1007 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "int_fast16_t", line: 1, size: 64, align: 64, baseType: !32)
!1008 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "int_fast32_t", line: 1, size: 64, align: 64, baseType: !32)
!1009 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "int_fast64_t", line: 1, size: 64, align: 64, baseType: !32)
!1010 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "uint_fast16_t", line: 1, size: 64, align: 64, baseType: !30)
!1011 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "uint_fast32_t", line: 1, size: 64, align: 64, baseType: !30)
!1012 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "uint_fast64_t", line: 1, size: 64, align: 64, baseType: !30)
!1013 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "intptr_t", line: 1, size: 64, align: 64, baseType: !32)
!1014 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "uintptr_t", line: 1, size: 64, align: 64, baseType: !30)
!1015 = !{ !1017, !1018, !1019, !1020, !1021, !1022, !1023, !1024, !1025, !1026, !1027, !1028, !1029, !1030, !1031, !1032, !1033, !1034, !1035, !1036, !1037, !1038, !1039, !1040 }
!1016 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "lconv", line: 1, size: 768, align: 64, elements: !1015, runtimeLang: DW_LANG_C_plus_plus)
!1017 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "decimal_point", line: 1, size: 64, align: 64, baseType: !44)
!1018 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "thousands_sep", line: 1, size: 64, align: 64, offset: 64, baseType: !44)
!1019 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "grouping", line: 1, size: 64, align: 64, offset: 128, baseType: !44)
!1020 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "int_curr_symbol", line: 1, size: 64, align: 64, offset: 192, baseType: !44)
!1021 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "currency_symbol", line: 1, size: 64, align: 64, offset: 256, baseType: !44)
!1022 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "mon_decimal_point", line: 1, size: 64, align: 64, offset: 320, baseType: !44)
!1023 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "mon_thousands_sep", line: 1, size: 64, align: 64, offset: 384, baseType: !44)
!1024 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "mon_grouping", line: 1, size: 64, align: 64, offset: 448, baseType: !44)
!1025 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "positive_sign", line: 1, size: 64, align: 64, offset: 512, baseType: !44)
!1026 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "negative_sign", line: 1, size: 64, align: 64, offset: 576, baseType: !44)
!1027 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "int_frac_digits", line: 1, size: 8, align: 8, offset: 640, baseType: !43)
!1028 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "frac_digits", line: 1, size: 8, align: 8, offset: 648, baseType: !43)
!1029 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "p_cs_precedes", line: 1, size: 8, align: 8, offset: 656, baseType: !43)
!1030 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "p_sep_by_space", line: 1, size: 8, align: 8, offset: 664, baseType: !43)
!1031 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "n_cs_precedes", line: 1, size: 8, align: 8, offset: 672, baseType: !43)
!1032 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "n_sep_by_space", line: 1, size: 8, align: 8, offset: 680, baseType: !43)
!1033 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "p_sign_posn", line: 1, size: 8, align: 8, offset: 688, baseType: !43)
!1034 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "n_sign_posn", line: 1, size: 8, align: 8, offset: 696, baseType: !43)
!1035 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "int_p_cs_precedes", line: 1, size: 8, align: 8, offset: 704, baseType: !43)
!1036 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "int_p_sep_by_space", line: 1, size: 8, align: 8, offset: 712, baseType: !43)
!1037 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "int_n_cs_precedes", line: 1, size: 8, align: 8, offset: 720, baseType: !43)
!1038 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "int_n_sep_by_space", line: 1, size: 8, align: 8, offset: 728, baseType: !43)
!1039 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "int_p_sign_posn", line: 1, size: 8, align: 8, offset: 736, baseType: !43)
!1040 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1016, name: "int_n_sign_posn", line: 1, size: 8, align: 8, offset: 744, baseType: !43)
!1041 = !{ !1043 }
!1042 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "sched_param", line: 1, size: 32, align: 32, elements: !1041, runtimeLang: DW_LANG_C_plus_plus)
!1043 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1042, name: "sched_priority", line: 1, size: 32, align: 32, baseType: !61)
!1044 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__cpu_mask", line: 1, size: 64, align: 64, baseType: !30)
!1045 = !{ !1047 }
!1046 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "cpu_set_t", line: 1, size: 1024, align: 64, elements: !1045, runtimeLang: DW_LANG_C_plus_plus)
!1047 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1046, name: "__bits", line: 1, size: 1024, align: 64, baseType: !329)
!1048 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "cpu_set_t", line: 1, size: 1024, align: 64, baseType: !1046)
!1049 = !{ !1051, !1052, !1053, !1054, !1055, !1056, !1057, !1058, !1059, !1060, !1061, !1062, !1063, !1064, !1065, !1066, !1067, !1068, !1069, !1070 }
!1050 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "timex", line: 1, size: 1664, align: 64, elements: !1049, runtimeLang: DW_LANG_C_plus_plus)
!1051 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "modes", line: 1, size: 32, align: 32, baseType: !97)
!1052 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "offset", line: 1, size: 64, align: 64, offset: 64, baseType: !32)
!1053 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "freq", line: 1, size: 64, align: 64, offset: 128, baseType: !32)
!1054 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "maxerror", line: 1, size: 64, align: 64, offset: 192, baseType: !32)
!1055 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "esterror", line: 1, size: 64, align: 64, offset: 256, baseType: !32)
!1056 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "status", line: 1, size: 32, align: 32, offset: 320, baseType: !61)
!1057 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "constant", line: 1, size: 64, align: 64, offset: 384, baseType: !32)
!1058 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "precision", line: 1, size: 64, align: 64, offset: 448, baseType: !32)
!1059 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "tolerance", line: 1, size: 64, align: 64, offset: 512, baseType: !32)
!1060 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "time", line: 1, size: 128, align: 64, offset: 576, baseType: !774)
!1061 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "tick", line: 1, size: 64, align: 64, offset: 704, baseType: !32)
!1062 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "ppsfreq", line: 1, size: 64, align: 64, offset: 768, baseType: !32)
!1063 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "jitter", line: 1, size: 64, align: 64, offset: 832, baseType: !32)
!1064 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "shift", line: 1, size: 32, align: 32, offset: 896, baseType: !61)
!1065 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "stabil", line: 1, size: 64, align: 64, offset: 960, baseType: !32)
!1066 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "jitcnt", line: 1, size: 64, align: 64, offset: 1024, baseType: !32)
!1067 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "calcnt", line: 1, size: 64, align: 64, offset: 1088, baseType: !32)
!1068 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "errcnt", line: 1, size: 64, align: 64, offset: 1152, baseType: !32)
!1069 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "stbcnt", line: 1, size: 64, align: 64, offset: 1216, baseType: !32)
!1070 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1050, name: "tai", line: 1, size: 32, align: 32, offset: 1280, baseType: !61)
!1071 = !{ !1073, !1074, !1075, !1076, !1077, !1078, !1079, !1080, !1081, !1082, !1083 }
!1072 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "tm", line: 1, size: 448, align: 64, elements: !1071, runtimeLang: DW_LANG_C_plus_plus)
!1073 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_sec", line: 1, size: 32, align: 32, baseType: !61)
!1074 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_min", line: 1, size: 32, align: 32, offset: 32, baseType: !61)
!1075 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_hour", line: 1, size: 32, align: 32, offset: 64, baseType: !61)
!1076 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_mday", line: 1, size: 32, align: 32, offset: 96, baseType: !61)
!1077 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_mon", line: 1, size: 32, align: 32, offset: 128, baseType: !61)
!1078 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_year", line: 1, size: 32, align: 32, offset: 160, baseType: !61)
!1079 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_wday", line: 1, size: 32, align: 32, offset: 192, baseType: !61)
!1080 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_yday", line: 1, size: 32, align: 32, offset: 224, baseType: !61)
!1081 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_isdst", line: 1, size: 32, align: 32, offset: 256, baseType: !61)
!1082 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_gmtoff", line: 1, size: 64, align: 64, offset: 320, baseType: !32)
!1083 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1072, name: "tm_zone", line: 1, size: 64, align: 64, offset: 384, baseType: !44)
!1084 = !{ !1086, !1087 }
!1085 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "itimerspec", line: 1, size: 256, align: 64, elements: !1084, runtimeLang: DW_LANG_C_plus_plus)
!1086 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1085, name: "it_interval", line: 1, size: 128, align: 64, baseType: !778)
!1087 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1085, name: "it_value", line: 1, size: 128, align: 64, offset: 128, baseType: !778)
!1088 = !{  }
!1089 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "sigevent", line: 1, align: 8, elements: !1088, runtimeLang: DW_LANG_C_plus_plus)
!1090 = !DICompositeType(tag: DW_TAG_array_type, size: 512, align: 64, baseType: !32, elements: !87)
!1091 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__jmp_buf", line: 1, size: 512, align: 64, baseType: !1090)
!1092 = !{ !1097, !1098, !1099, !1101 }
!1093 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_pthread_cleanup_buffer", line: 1, size: 256, align: 64, elements: !1092, runtimeLang: DW_LANG_C_plus_plus)
!1094 = !{ null, !17 }
!1095 = !DISubroutineType(types: !1094)
!1096 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1095)
!1097 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1093, name: "__routine", line: 1, size: 64, align: 64, baseType: !1096)
!1098 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1093, name: "__arg", line: 1, size: 64, align: 64, offset: 64, baseType: !17)
!1099 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1093, name: "__canceltype", line: 1, size: 32, align: 32, offset: 128, baseType: !61)
!1100 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1093)
!1101 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1093, name: "__prev", line: 1, size: 64, align: 64, offset: 192, baseType: !1100)
!1102 = !{ !1109, !1111 }
!1103 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pthread_unwind_buf_t", line: 1, size: 832, align: 64, elements: !1102, runtimeLang: DW_LANG_C_plus_plus)
!1104 = !{ !1106, !1107 }
!1105 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !1103, name: "_ZN22__pthread_unwind_buf_tUt_E", line: 1, size: 576, align: 64, elements: !1104, runtimeLang: DW_LANG_C_plus_plus)
!1106 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1105, name: "__cancel_jmp_buf", line: 1, size: 512, align: 64, baseType: !1090)
!1107 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1105, name: "__mask_was_saved", line: 1, size: 32, align: 32, offset: 512, baseType: !61)
!1108 = !DICompositeType(tag: DW_TAG_array_type, size: 576, align: 64, baseType: !1105, elements: !377)
!1109 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1103, name: "__cancel_jmp_buf", line: 1, size: 576, align: 64, baseType: !1108)
!1110 = !DICompositeType(tag: DW_TAG_array_type, size: 256, align: 64, baseType: !17, elements: !69)
!1111 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1103, name: "__pad", line: 1, size: 256, align: 64, offset: 576, baseType: !1110)
!1112 = !{ !1114, !1115, !1116, !1117 }
!1113 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__pthread_cleanup_frame", line: 1, size: 192, align: 64, elements: !1112, runtimeLang: DW_LANG_C_plus_plus)
!1114 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1113, name: "__cancel_routine", line: 1, size: 64, align: 64, baseType: !1096)
!1115 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1113, name: "__cancel_arg", line: 1, size: 64, align: 64, offset: 64, baseType: !17)
!1116 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1113, name: "__do_it", line: 1, size: 32, align: 32, offset: 128, baseType: !61)
!1117 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1113, name: "__cancel_type", line: 1, size: 32, align: 32, offset: 160, baseType: !61)
!1118 = !{ !1120, !1121, !1122, !1123 }
!1119 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "__pthread_cleanup_class", line: 1, size: 192, align: 64, elements: !1118, runtimeLang: DW_LANG_C_plus_plus)
!1120 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1119, name: "__cancel_routine", line: 1, size: 64, align: 64, baseType: !1096)
!1121 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1119, name: "__cancel_arg", line: 1, size: 64, align: 64, offset: 64, baseType: !17)
!1122 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1119, name: "__do_it", line: 1, size: 32, align: 32, offset: 128, baseType: !61)
!1123 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1119, name: "__cancel_type", line: 1, size: 32, align: 32, offset: 160, baseType: !61)
!1124 = !{  }
!1125 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "__jmp_buf_tag", line: 1, align: 8, elements: !1124, runtimeLang: DW_LANG_C_plus_plus)
!1126 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gthread_t", line: 1, size: 64, align: 64, baseType: !30)
!1127 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gthread_key_t", line: 1, size: 32, align: 32, baseType: !97)
!1128 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gthread_once_t", line: 1, size: 32, align: 32, baseType: !61)
!1129 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gthread_mutex_t", line: 1, size: 320, align: 64, baseType: !873)
!1130 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gthread_recursive_mutex_t", line: 1, size: 320, align: 64, baseType: !873)
!1131 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gthread_cond_t", line: 1, size: 384, align: 64, baseType: !882)
!1132 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__gthread_time_t", line: 1, size: 128, align: 64, baseType: !778)
!1133 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_Atomic_word", line: 1, size: 32, align: 32, baseType: !61)
!1134 = !{ !1136, !1137 }
!1135 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_G_fpos_t", line: 1, size: 128, align: 64, elements: !1134, runtimeLang: DW_LANG_C_plus_plus)
!1136 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1135, name: "__pos", line: 1, size: 64, align: 64, baseType: !32)
!1137 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1135, name: "__state", line: 1, size: 64, align: 32, offset: 64, baseType: !945)
!1138 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fpos_t", line: 1, size: 128, align: 64, baseType: !1135)
!1139 = !{ !1141, !1142 }
!1140 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_G_fpos64_t", line: 1, size: 128, align: 64, elements: !1139, runtimeLang: DW_LANG_C_plus_plus)
!1141 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1140, name: "__pos", line: 1, size: 64, align: 64, baseType: !32)
!1142 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1140, name: "__state", line: 1, size: 64, align: 32, offset: 64, baseType: !945)
!1143 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__fpos64_t", line: 1, size: 128, align: 64, baseType: !1140)
!1144 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_IO_lock_t", line: 1, align: 1, baseType: !16)
!1145 = !{ !32, !17, !44, !30 }
!1146 = !DISubroutineType(types: !1145)
!1147 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "cookie_read_function_t", line: 1, size: 8, align: 1, baseType: !1146)
!1148 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__ssize_t", line: 1, size: 64, align: 64, baseType: !32)
!1149 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "size_t", line: 1, size: 64, align: 64, baseType: !30)
!1150 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "cookie_write_function_t", line: 1, size: 8, align: 1, baseType: !1146)
!1151 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "__off64_t", line: 1, size: 64, align: 64, baseType: !32)
!1152 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !32)
!1153 = !{ !61, !17, !1152, !61 }
!1154 = !DISubroutineType(types: !1153)
!1155 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "cookie_seek_function_t", line: 1, size: 8, align: 1, baseType: !1154)
!1156 = !{ !61, !17 }
!1157 = !DISubroutineType(types: !1156)
!1158 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "cookie_close_function_t", line: 1, size: 8, align: 1, baseType: !1157)
!1159 = !{ !1162, !1163, !1165, !1167 }
!1160 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_IO_cookie_io_functions_t", line: 1, size: 256, align: 64, elements: !1159, runtimeLang: DW_LANG_C_plus_plus)
!1161 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1146)
!1162 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1160, name: "read", line: 1, size: 64, align: 64, baseType: !1161)
!1163 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1160, name: "write", line: 1, size: 64, align: 64, offset: 64, baseType: !1161)
!1164 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1154)
!1165 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1160, name: "seek", line: 1, size: 64, align: 64, offset: 128, baseType: !1164)
!1166 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1157)
!1167 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1160, name: "close", line: 1, size: 64, align: 64, offset: 192, baseType: !1166)
!1168 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "cookie_io_functions_t", line: 1, size: 256, align: 64, baseType: !1160)
!1169 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "fpos_t", line: 1, size: 128, align: 64, baseType: !1135)
!1170 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "fpos64_t", line: 1, size: 128, align: 64, baseType: !1140)
!1171 = !{  }
!1172 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "obstack", line: 1, align: 8, elements: !1171, runtimeLang: DW_LANG_C_plus_plus)
!1173 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "error_t", line: 1, size: 32, align: 32, baseType: !61)
!1174 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "wctype_t", line: 1, size: 64, align: 64, baseType: !30)
!1175 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "wctrans_t", line: 1, size: 64, align: 64, baseType: !62)
!1176 = !{  }
!1177 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T1DFunction", line: 1, size: 64, align: 64, elements: !1176, runtimeLang: DW_LANG_C_plus_plus)
!1178 = !{  }
!1179 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T2DFunction", line: 1, size: 64, align: 64, elements: !1178, runtimeLang: DW_LANG_C_plus_plus)
!1180 = !{  }
!1181 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T3DFunction", line: 1, size: 64, align: 64, elements: !1180, runtimeLang: DW_LANG_C_plus_plus)
!1182 = !{ !1184, !1186, !1187, !1191 }
!1183 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T2D_fix1", line: 1, size: 192, align: 64, elements: !1182, runtimeLang: DW_LANG_C_plus_plus)
!1184 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1183, name: "T1DFunction", line: 1, size: 64, align: 64, baseType: !1177)
!1185 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1179)
!1186 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1183, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1185)
!1187 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1183, name: "x", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1188 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1183)
!1189 = !{ !731, !1188, !731 }
!1190 = !DISubroutineType(types: !1189)
!1191 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1183, name: "call", line: 41, size: 8, align: 1, baseType: !1190)
!1192 = !{ !1194, !1195, !1196 }
!1193 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T2D_fix2", line: 1, size: 192, align: 64, elements: !1192, runtimeLang: DW_LANG_C_plus_plus)
!1194 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1193, name: "T1DFunction", line: 1, size: 64, align: 64, baseType: !1177)
!1195 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1193, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1185)
!1196 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1193, name: "y", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1197 = !{ !1199, !1201, !1202 }
!1198 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T3D_fix1", line: 1, size: 192, align: 64, elements: !1197, runtimeLang: DW_LANG_C_plus_plus)
!1199 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1198, name: "T2DFunction", line: 1, size: 64, align: 64, baseType: !1179)
!1200 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1181)
!1201 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1198, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1200)
!1202 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1198, name: "x", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1203 = !{ !1205, !1206, !1207 }
!1204 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T3D_fix2", line: 1, size: 192, align: 64, elements: !1203, runtimeLang: DW_LANG_C_plus_plus)
!1205 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1204, name: "T2DFunction", line: 1, size: 64, align: 64, baseType: !1179)
!1206 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1204, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1200)
!1207 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1204, name: "y", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1208 = !{ !1210, !1211, !1212, !1216 }
!1209 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T3D_fix3", line: 1, size: 192, align: 64, elements: !1208, runtimeLang: DW_LANG_C_plus_plus)
!1210 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1209, name: "T2DFunction", line: 1, size: 64, align: 64, baseType: !1179)
!1211 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1209, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1200)
!1212 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1209, name: "z", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1213 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1209)
!1214 = !{ !731, !1213, !731, !731 }
!1215 = !DISubroutineType(types: !1214)
!1216 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1209, name: "call", line: 83, size: 8, align: 1, baseType: !1215)
!1217 = !{ !1219, !1220, !1221, !1222 }
!1218 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T3D_fix12", line: 1, size: 256, align: 64, elements: !1217, runtimeLang: DW_LANG_C_plus_plus)
!1219 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1218, name: "T1DFunction", line: 1, size: 64, align: 64, baseType: !1177)
!1220 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1218, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1200)
!1221 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1218, name: "x", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1222 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1218, name: "y", line: 1, size: 64, align: 64, offset: 192, baseType: !731)
!1223 = !{ !1225, !1226, !1227, !1228 }
!1224 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T3D_fix13", line: 1, size: 256, align: 64, elements: !1223, runtimeLang: DW_LANG_C_plus_plus)
!1225 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1224, name: "T1DFunction", line: 1, size: 64, align: 64, baseType: !1177)
!1226 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1224, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1200)
!1227 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1224, name: "x", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1228 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1224, name: "z", line: 1, size: 64, align: 64, offset: 192, baseType: !731)
!1229 = !{ !1231, !1232, !1233, !1234 }
!1230 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "T3D_fix23", line: 1, size: 256, align: 64, elements: !1229, runtimeLang: DW_LANG_C_plus_plus)
!1231 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1230, name: "T1DFunction", line: 1, size: 64, align: 64, baseType: !1177)
!1232 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1230, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1200)
!1233 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1230, name: "y", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1234 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1230, name: "z", line: 1, size: 64, align: 64, offset: 192, baseType: !731)
!1235 = !{ !1237, !1238, !1239, !1240, !1241, !1245 }
!1236 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "Tinty_f2D", line: 1, size: 320, align: 64, elements: !1235, runtimeLang: DW_LANG_C_plus_plus)
!1237 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1236, name: "T1DFunction", line: 1, size: 64, align: 64, baseType: !1177)
!1238 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1236, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1185)
!1239 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1236, name: "ymin2D", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1240 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1236, name: "ymax2D", line: 1, size: 64, align: 64, offset: 192, baseType: !731)
!1241 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1236, name: "the_absacc2D", line: 1, size: 64, align: 64, offset: 256, baseType: !731)
!1242 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1236)
!1243 = !{ !731, !1242, !731 }
!1244 = !DISubroutineType(types: !1243)
!1245 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1236, name: "call", line: 226, size: 8, align: 1, baseType: !1244)
!1246 = !{ !1248, !1249, !1250, !1251, !1252, !1253, !1254, !1258 }
!1247 = !DICompositeType(tag: DW_TAG_class_type, file: !3, scope: !10, name: "Tintxy_f3D", line: 1, size: 448, align: 64, elements: !1246, runtimeLang: DW_LANG_C_plus_plus)
!1248 = !DIDerivedType(tag: DW_TAG_inheritance, file: !3, scope: !1247, name: "T1DFunction", line: 1, size: 64, align: 64, baseType: !1177)
!1249 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1247, name: "f", line: 1, size: 64, align: 64, offset: 64, baseType: !1200)
!1250 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1247, name: "xmin3D", line: 1, size: 64, align: 64, offset: 128, baseType: !731)
!1251 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1247, name: "xmax3D", line: 1, size: 64, align: 64, offset: 192, baseType: !731)
!1252 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1247, name: "ymin3D", line: 1, size: 64, align: 64, offset: 256, baseType: !731)
!1253 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1247, name: "ymax3D", line: 1, size: 64, align: 64, offset: 320, baseType: !731)
!1254 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1247, name: "the_absacc3D", line: 1, size: 64, align: 64, offset: 384, baseType: !731)
!1255 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1247)
!1256 = !{ !731, !1255, !731 }
!1257 = !DISubroutineType(types: !1256)
!1258 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1247, name: "call", line: 244, size: 8, align: 1, baseType: !1257)
!1259 = !{ !1261 }
!1260 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1259, runtimeLang: DW_LANG_C_plus_plus)
!1261 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1260, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1262 = !{  }
!1263 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1262, runtimeLang: DW_LANG_C_plus_plus)
!1264 = !{ !1266 }
!1265 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1264, runtimeLang: DW_LANG_C_plus_plus)
!1266 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1265, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1267 = !{  }
!1268 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1267, runtimeLang: DW_LANG_C_plus_plus)
!1269 = !{ !1271 }
!1270 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1269, runtimeLang: DW_LANG_C_plus_plus)
!1271 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1270, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1272 = !{  }
!1273 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1272, runtimeLang: DW_LANG_C_plus_plus)
!1274 = !{ !1276 }
!1275 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1274, runtimeLang: DW_LANG_C_plus_plus)
!1276 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1275, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1277 = !{  }
!1278 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1277, runtimeLang: DW_LANG_C_plus_plus)
!1279 = !{ !1281 }
!1280 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1279, runtimeLang: DW_LANG_C_plus_plus)
!1281 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1280, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1282 = !{  }
!1283 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1282, runtimeLang: DW_LANG_C_plus_plus)
!1284 = !{ !1286 }
!1285 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1284, runtimeLang: DW_LANG_C_plus_plus)
!1286 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1285, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1287 = !{  }
!1288 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1287, runtimeLang: DW_LANG_C_plus_plus)
!1289 = !{ !1291 }
!1290 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1289, runtimeLang: DW_LANG_C_plus_plus)
!1291 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1290, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1292 = !{  }
!1293 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1292, runtimeLang: DW_LANG_C_plus_plus)
!1294 = !{ !1296 }
!1295 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1294, runtimeLang: DW_LANG_C_plus_plus)
!1296 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1295, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1297 = !{  }
!1298 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1297, runtimeLang: DW_LANG_C_plus_plus)
!1299 = !{ !1301 }
!1300 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1299, runtimeLang: DW_LANG_C_plus_plus)
!1301 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1300, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1302 = !{  }
!1303 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1302, runtimeLang: DW_LANG_C_plus_plus)
!1304 = !{ !1306 }
!1305 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1304, runtimeLang: DW_LANG_C_plus_plus)
!1306 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1305, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1307 = !{  }
!1308 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1307, runtimeLang: DW_LANG_C_plus_plus)
!1309 = !{ !1311 }
!1310 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1309, runtimeLang: DW_LANG_C_plus_plus)
!1311 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1310, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1312 = !{  }
!1313 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1312, runtimeLang: DW_LANG_C_plus_plus)
!1314 = !{ !1316 }
!1315 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1314, runtimeLang: DW_LANG_C_plus_plus)
!1316 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1315, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1317 = !{  }
!1318 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1317, runtimeLang: DW_LANG_C_plus_plus)
!1319 = !{ !1321 }
!1320 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1319, runtimeLang: DW_LANG_C_plus_plus)
!1321 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1320, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1322 = !{  }
!1323 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1322, runtimeLang: DW_LANG_C_plus_plus)
!1324 = !{ !1326 }
!1325 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1324, runtimeLang: DW_LANG_C_plus_plus)
!1326 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1325, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1327 = !{  }
!1328 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1327, runtimeLang: DW_LANG_C_plus_plus)
!1329 = !{ !1331 }
!1330 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1329, runtimeLang: DW_LANG_C_plus_plus)
!1331 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1330, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1332 = !{  }
!1333 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1332, runtimeLang: DW_LANG_C_plus_plus)
!1334 = !{ !1336 }
!1335 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Save_errno", line: 1, size: 32, align: 32, elements: !1334, runtimeLang: DW_LANG_C_plus_plus)
!1336 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1335, name: "_M_errno", line: 1, size: 32, align: 32, baseType: !61)
!1337 = !{  }
!1338 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Range_chk", line: 1, size: 8, align: 8, elements: !1337, runtimeLang: DW_LANG_C_plus_plus)
!1339 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_Tag", line: 1, size: 8, align: 8, baseType: !440)
!1340 = !DIDerivedType(tag: DW_TAG_typedef, file: !3, scope: !10, name: "_Integral", line: 1, size: 8, align: 8, baseType: !38)
!1341 = !{  }
!1342 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Iter_less_iter", line: 1, size: 8, align: 8, elements: !1341, runtimeLang: DW_LANG_C_plus_plus)
!1343 = !{  }
!1344 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Iter_less_val", line: 1, size: 8, align: 8, elements: !1343, runtimeLang: DW_LANG_C_plus_plus)
!1345 = !{  }
!1346 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Val_less_iter", line: 1, size: 8, align: 8, elements: !1345, runtimeLang: DW_LANG_C_plus_plus)
!1347 = !{  }
!1348 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Iter_equal_to_iter", line: 1, size: 8, align: 8, elements: !1347, runtimeLang: DW_LANG_C_plus_plus)
!1349 = !{  }
!1350 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, scope: !10, name: "_Iter_equal_to_val", line: 1, size: 8, align: 8, elements: !1349, runtimeLang: DW_LANG_C_plus_plus)
!1351 = !{ !1353 }
!1352 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "T1DFunction", size: 64, align: 64, elements: !1351)
!1353 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1352, name: "__vptr", size: 64, align: 64, baseType: !126)
!1354 = !DIDerivedType(tag: DW_TAG_reference_type, size: 64, align: 64, baseType: !1352)
!1355 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !731)
!1356 = !{ null, !1354, !731, !731, !1355, !62, !61 }
!1357 = !DISubroutineType(types: !1356)
!1358 = !{ !1367, !1369, !1371, !1373, !1375, !1377 }
!1359 = distinct !DISubprogram(file: !3, scope: !10, name: "trapez", line: 46, type: !1357, spFlags: 12, unit: !10, scopeLine: 46)
!1360 = !DILocation(scope: !1359)
!1361 = !DILexicalBlock(file: !3, scope: !1359, line: 46, column: 1)
!1362 = !DILocation(scope: !1361)
!1363 = !DILexicalBlock(file: !3, scope: !1361, line: 52, column: 1)
!1364 = !DILocation(scope: !1363)
!1365 = !DILocalVariable(scope: !1361, name: "func", file: !3, type: !1354)
!1366 = !DIExpression()
!1367 = !DILocalVariable(scope: !1359, name: "func", arg: 1, file: !3, type: !1354)
!1368 = !DILocalVariable(scope: !1361, name: "a", file: !3, type: !731)
!1369 = !DILocalVariable(scope: !1359, name: "a", arg: 2, file: !3, type: !731)
!1370 = !DILocalVariable(scope: !1361, name: "b", file: !3, type: !731)
!1371 = !DILocalVariable(scope: !1359, name: "b", arg: 3, file: !3, type: !731)
!1372 = !DILocalVariable(scope: !1361, name: "S", file: !3, type: !1355)
!1373 = !DILocalVariable(scope: !1359, name: "S", arg: 4, file: !3, type: !1355)
!1374 = !DILocalVariable(scope: !1361, name: "it", file: !3, type: !62)
!1375 = !DILocalVariable(scope: !1359, name: "it", arg: 5, file: !3, type: !62)
!1376 = !DILocalVariable(scope: !1361, name: "n", file: !3, type: !61)
!1377 = !DILocalVariable(scope: !1359, name: "n", arg: 6, file: !3, type: !61)
!1378 = !DILocation(line: 48, column: 1, scope: !1361)
!1379 = !DILocation(line: 49, column: 1, scope: !1361)
!1380 = !DILocation(line: 50, column: 1, scope: !1361)
!1381 = !DILocation(line: 52, column: 1, scope: !1363)
!1382 = !DILocation(line: 53, column: 1, scope: !1363)
!1383 = !DILocalVariable(scope: !1363, name: "delta", file: !3, type: !731)
!1384 = !DILocation(line: 54, column: 1, scope: !1363)
!1385 = !DILocalVariable(scope: !1363, name: "x", file: !3, type: !731)
!1386 = !DILocation(line: 55, column: 1, scope: !1363)
!1387 = !DILocalVariable(scope: !1363, name: "sum", file: !3, type: !731)
!1388 = !DILocation(line: 56, column: 1, scope: !1363)
!1389 = !DILocalVariable(scope: !1361, name: "j", file: !3, type: !61)
!1390 = !DILocation(line: 57, column: 1, scope: !1363)
!1391 = !DILocation(line: 58, column: 1, scope: !1363)
!1392 = !DILocation(line: 59, column: 1, scope: !1363)
!1393 = !DILocation(line: 60, column: 1, scope: !1363)
!1394 = !DILocation(line: 61, column: 1, scope: !1363)
!1395 = !DILocation(line: 64, column: 1, scope: !1361)
!1396 = !{ !"PGI C[++] TBAA" }
!1397 = !{ !"omnipotent char", !1396, i64 0 }
!1398 = !{ !"<T>*", !1397, i64 0 }
!1399 = !{ !1397, !1397, i64 0 }
!1400 = !{ !"int", !1397, i64 0 }
!1401 = !{ !1400, !1400, i64 0 }
!1402 = !{ !"double", !1397, i64 0 }
!1403 = !{ !1402, !1402, i64 0 }
!1404 = !DIDerivedType(tag: DW_TAG_reference_type, size: 64, align: 64, baseType: !731)
!1405 = !{ null, !1355, !1355, !61, !731, !1404, !1404 }
!1406 = !DISubroutineType(types: !1405)
!1407 = !{ !1419, !1421, !1423, !1425, !1427, !1429 }
!1408 = distinct !DISubprogram(file: !3, scope: !10, name: "polint", line: 70, type: !1406, spFlags: 12, unit: !10, scopeLine: 70)
!1409 = !DILocation(scope: !1408)
!1410 = !DILexicalBlock(file: !3, scope: !1408, line: 70, column: 1)
!1411 = !DILocation(scope: !1410)
!1412 = !DILexicalBlock(file: !3, scope: !1410, line: 76, column: 1)
!1413 = !DILocation(scope: !1412)
!1414 = !DILexicalBlock(file: !3, scope: !1410, line: 87, column: 1)
!1415 = !DILocation(scope: !1414)
!1416 = !DILexicalBlock(file: !3, scope: !1414, line: 88, column: 1)
!1417 = !DILocation(scope: !1416)
!1418 = !DILocalVariable(scope: !1410, name: "xa", file: !3, type: !1355)
!1419 = !DILocalVariable(scope: !1408, name: "xa", arg: 1, file: !3, type: !1355)
!1420 = !DILocalVariable(scope: !1410, name: "ya", file: !3, type: !1355)
!1421 = !DILocalVariable(scope: !1408, name: "ya", arg: 2, file: !3, type: !1355)
!1422 = !DILocalVariable(scope: !1410, name: "n", file: !3, type: !61)
!1423 = !DILocalVariable(scope: !1408, name: "n", arg: 3, file: !3, type: !61)
!1424 = !DILocalVariable(scope: !1410, name: "x", file: !3, type: !731)
!1425 = !DILocalVariable(scope: !1408, name: "x", arg: 4, file: !3, type: !731)
!1426 = !DILocalVariable(scope: !1410, name: "y", file: !3, type: !1404)
!1427 = !DILocalVariable(scope: !1408, name: "y", arg: 5, file: !3, type: !1404)
!1428 = !DILocalVariable(scope: !1410, name: "dy", file: !3, type: !1404)
!1429 = !DILocalVariable(scope: !1408, name: "dy", arg: 6, file: !3, type: !1404)
!1430 = !DILocation(line: 74, column: 1, scope: !1410)
!1431 = !DILocalVariable(scope: !1410, name: "ns", file: !3, type: !61)
!1432 = !DILocation(line: 75, column: 1, scope: !1410)
!1433 = !DILocalVariable(scope: !1410, name: "dif", file: !3, type: !731)
!1434 = !DILocation(line: 76, column: 1, scope: !1412)
!1435 = !DILocalVariable(scope: !1412, name: "dift", file: !3, type: !731)
!1436 = !DICompositeType(tag: DW_TAG_array_type, size: 1280, align: 64, baseType: !731, elements: !995)
!1437 = !DILocalVariable(scope: !1410, name: "D", file: !3, type: !1436)
!1438 = !DILocalVariable(scope: !1410, name: "C", file: !3, type: !1436)
!1439 = !DILocation(line: 85, column: 1, scope: !1410)
!1440 = !DILocation(line: 86, column: 1, scope: !1410)
!1441 = !DILocation(line: 106, column: 1, scope: !1410)
!1442 = !DILocation(line: 87, column: 1, scope: !1414)
!1443 = !DILocation(line: 88, column: 1, scope: !1416)
!1444 = !DILocalVariable(scope: !1416, name: "h0", file: !3, type: !731)
!1445 = !DILocalVariable(scope: !1416, name: "hp", file: !3, type: !731)
!1446 = !DILocalVariable(scope: !1416, name: "W", file: !3, type: !731)
!1447 = !DILocalVariable(scope: !1416, name: "den", file: !3, type: !731)
!1448 = !{ !"long", !1397, i64 0 }
!1449 = !{ !1448, !1448, i64 0 }
!1450 = !{ !"long long", !1397, i64 0 }
!1451 = !{ !1450, !1450, i64 0 }
!1452 = !{ !1478, !1480, !1482, !1484, !1486, !1488 }
!1453 = distinct !DISubprogram(file: !3, scope: !10, name: "ratint", line: 112, type: !1406, spFlags: 12, unit: !10, scopeLine: 112)
!1454 = !DILocation(scope: !1453)
!1455 = !DILexicalBlock(file: !3, scope: !1453, line: 112, column: 1)
!1456 = !DILocation(scope: !1455)
!1457 = !DILexicalBlock(file: !3, scope: !1455, line: 1, column: 1)
!1458 = !DILocation(scope: !1457)
!1459 = !DILexicalBlock(file: !3, scope: !1457, line: 1, column: 1)
!1460 = !DILocation(scope: !1459)
!1461 = !DILexicalBlock(file: !3, scope: !1459, line: 1, column: 1)
!1462 = !DILocation(scope: !1461)
!1463 = !DILexicalBlock(file: !3, scope: !1461, line: 1, column: 1)
!1464 = !DILocation(scope: !1463)
!1465 = !DILexicalBlock(file: !3, scope: !1457, line: 1, column: 1)
!1466 = !DILocation(scope: !1465)
!1467 = !DILexicalBlock(file: !3, scope: !1455, line: 1, column: 1)
!1468 = !DILocation(scope: !1467)
!1469 = !DILexicalBlock(file: !3, scope: !1467, line: 1, column: 1)
!1470 = !DILocation(scope: !1469)
!1471 = !DILexicalBlock(file: !3, scope: !1469, line: 1, column: 1)
!1472 = !DILocation(scope: !1471)
!1473 = !DILexicalBlock(file: !3, scope: !1471, line: 1, column: 1)
!1474 = !DILocation(scope: !1473)
!1475 = !DILexicalBlock(file: !3, scope: !1467, line: 1, column: 1)
!1476 = !DILocation(scope: !1475)
!1477 = !DILocalVariable(scope: !1455, name: "xa", file: !3, type: !1355)
!1478 = !DILocalVariable(scope: !1453, name: "xa", arg: 1, file: !3, type: !1355)
!1479 = !DILocalVariable(scope: !1455, name: "ya", file: !3, type: !1355)
!1480 = !DILocalVariable(scope: !1453, name: "ya", arg: 2, file: !3, type: !1355)
!1481 = !DILocalVariable(scope: !1455, name: "n", file: !3, type: !61)
!1482 = !DILocalVariable(scope: !1453, name: "n", arg: 3, file: !3, type: !61)
!1483 = !DILocalVariable(scope: !1455, name: "x", file: !3, type: !731)
!1484 = !DILocalVariable(scope: !1453, name: "x", arg: 4, file: !3, type: !731)
!1485 = !DILocalVariable(scope: !1455, name: "y", file: !3, type: !1404)
!1486 = !DILocalVariable(scope: !1453, name: "y", arg: 5, file: !3, type: !1404)
!1487 = !DILocalVariable(scope: !1455, name: "dy", file: !3, type: !1404)
!1488 = !DILocalVariable(scope: !1453, name: "dy", arg: 6, file: !3, type: !1404)
!1489 = !DILocation(line: 113, column: 1, scope: !1455)
!1490 = !DILocalVariable(scope: !1455, name: "ns", file: !3, type: !61)
!1491 = !DILocation(line: 118, column: 1, scope: !1455)
!1492 = !DILocalVariable(scope: !1455, name: "hh", file: !3, type: !731)
!1493 = !DILocation(line: 119, column: 1, scope: !1455)
!1494 = !DILocalVariable(scope: !1455, name: "i", file: !3, type: !61)
!1495 = !DILocation(line: 120, column: 1, scope: !1455)
!1496 = !DILocalVariable(scope: !1455, name: "h", file: !3, type: !731)
!1497 = !DILocation(line: 121, column: 1, scope: !1455)
!1498 = !DILocation(line: 122, column: 1, scope: !1455)
!1499 = !DILocation(line: 123, column: 1, scope: !1455)
!1500 = !DILocation(line: 124, column: 1, scope: !1455)
!1501 = !DILocation(line: 125, column: 1, scope: !1455)
!1502 = !DILocation(line: 126, column: 1, scope: !1455)
!1503 = !DILocation(line: 127, column: 1, scope: !1455)
!1504 = !DILocation(line: 129, column: 1, scope: !1455)
!1505 = !DILocalVariable(scope: !1455, name: "C", file: !3, type: !1436)
!1506 = !DILocation(line: 130, column: 1, scope: !1455)
!1507 = !DILocalVariable(scope: !1455, name: "D", file: !3, type: !1436)
!1508 = !DILocation(line: 131, column: 1, scope: !1455)
!1509 = !DILocation(line: 133, column: 1, scope: !1455)
!1510 = !DILocation(line: 135, column: 1, scope: !1455)
!1511 = !DILocalVariable(scope: !1455, name: "m1", file: !3, type: !61)
!1512 = !DILocation(line: 136, column: 1, scope: !1455)
!1513 = !DILocation(line: 137, column: 1, scope: !1455)
!1514 = !DILocalVariable(scope: !1455, name: "w", file: !3, type: !731)
!1515 = !DILocation(line: 138, column: 1, scope: !1455)
!1516 = !DILocation(line: 139, column: 1, scope: !1455)
!1517 = !DILocalVariable(scope: !1455, name: "t", file: !3, type: !731)
!1518 = !DILocation(line: 140, column: 1, scope: !1455)
!1519 = !DILocalVariable(scope: !1455, name: "dd", file: !3, type: !731)
!1520 = !DILocation(line: 141, column: 1, scope: !1455)
!1521 = !DILocation(line: 142, column: 1, scope: !1455)
!1522 = !DILocation(line: 568, column: 1, scope: !1455)
!1523 = distinct !DIGlobalVariable(scope: !11, name: "_ZSt4cerr", file: !3, type: !285)
!1524 = !DIGlobalVariableExpression(var: !1523, expr: !1366)
!1525 = !DILocation(line: 570, column: 1, scope: !1455)
!1526 = !DILocation(line: 143, column: 1, scope: !1455)
!1527 = !DILocation(line: 146, column: 1, scope: !1455)
!1528 = !DILocation(line: 147, column: 1, scope: !1455)
!1529 = !DILocation(line: 148, column: 1, scope: !1455)
!1530 = !DILocation(line: 149, column: 1, scope: !1455)
!1531 = !DILocation(line: 151, column: 1, scope: !1455)
!1532 = !DILocation(line: 152, column: 1, scope: !1455)
!1533 = !DILocation(line: 153, column: 1, scope: !1455)
!1534 = !{ !731, !1354, !731, !731, !731 }
!1535 = !DISubroutineType(types: !1534)
!1536 = !{ !1566, !1568, !1570, !1572 }
!1537 = distinct !DISubprogram(file: !3, scope: !10, name: "Romberg_simple", line: 156, type: !1535, spFlags: 8, unit: !10, scopeLine: 156)
!1538 = !DILocation(scope: !1537)
!1539 = !DILexicalBlock(file: !3, scope: !1537, line: 156, column: 1)
!1540 = !DILocation(scope: !1539)
!1541 = !DILexicalBlock(file: !3, scope: !1539, line: 1, column: 1)
!1542 = !DILocation(scope: !1541)
!1543 = !DILexicalBlock(file: !3, scope: !1541, line: 52, column: 1)
!1544 = !DILocation(scope: !1543)
!1545 = !DILexicalBlock(file: !3, scope: !1539, line: 1, column: 1)
!1546 = !DILocation(scope: !1545)
!1547 = !DILexicalBlock(file: !3, scope: !1545, line: 76, column: 1)
!1548 = !DILocation(scope: !1547)
!1549 = !DILexicalBlock(file: !3, scope: !1545, line: 87, column: 1)
!1550 = !DILocation(scope: !1549)
!1551 = !DILexicalBlock(file: !3, scope: !1549, line: 88, column: 1)
!1552 = !DILocation(scope: !1551)
!1553 = !DILexicalBlock(file: !3, scope: !1539, line: 1, column: 1)
!1554 = !DILocation(scope: !1553)
!1555 = !DILexicalBlock(file: !3, scope: !1553, line: 52, column: 1)
!1556 = !DILocation(scope: !1555)
!1557 = !DILexicalBlock(file: !3, scope: !1539, line: 1, column: 1)
!1558 = !DILocation(scope: !1557)
!1559 = !DILexicalBlock(file: !3, scope: !1557, line: 76, column: 1)
!1560 = !DILocation(scope: !1559)
!1561 = !DILexicalBlock(file: !3, scope: !1557, line: 87, column: 1)
!1562 = !DILocation(scope: !1561)
!1563 = !DILexicalBlock(file: !3, scope: !1561, line: 88, column: 1)
!1564 = !DILocation(scope: !1563)
!1565 = !DILocalVariable(scope: !1539, name: "func", file: !3, type: !1354)
!1566 = !DILocalVariable(scope: !1537, name: "func", arg: 1, file: !3, type: !1354)
!1567 = !DILocalVariable(scope: !1539, name: "a", file: !3, type: !731)
!1568 = !DILocalVariable(scope: !1537, name: "a", arg: 2, file: !3, type: !731)
!1569 = !DILocalVariable(scope: !1539, name: "b", file: !3, type: !731)
!1570 = !DILocalVariable(scope: !1537, name: "b", arg: 3, file: !3, type: !731)
!1571 = !DILocalVariable(scope: !1539, name: "absacc", file: !3, type: !731)
!1572 = !DILocalVariable(scope: !1537, name: "absacc", arg: 4, file: !3, type: !731)
!1573 = !DILocation(line: 157, column: 1, scope: !1539)
!1574 = !DILocalVariable(scope: !1539, name: "it", file: !3, type: !61)
!1575 = !DILocation(line: 162, column: 1, scope: !1539)
!1576 = !DISubrange(count: 9)
!1577 = !{ !1576 }
!1578 = !DICompositeType(tag: DW_TAG_array_type, size: 576, align: 64, baseType: !731, elements: !1577)
!1579 = !DILocalVariable(scope: !1539, name: "H", file: !3, type: !1578)
!1580 = !DILocation(line: 163, column: 1, scope: !1539)
!1581 = !DILocalVariable(scope: !1539, name: "dresult", file: !3, type: !731)
!1582 = !DILocation(line: 164, column: 1, scope: !1539)
!1583 = !DILocalVariable(scope: !1539, name: "j", file: !3, type: !61)
!1584 = !DILocalVariable(scope: !1539, name: "S", file: !3, type: !1578)
!1585 = !DILocation(line: 165, column: 1, scope: !1539)
!1586 = !DILocation(line: 49, column: 1, scope: !1539)
!1587 = !DILocation(line: 50, column: 1, scope: !1539)
!1588 = !DILocation(line: 52, column: 1, scope: !1555)
!1589 = !DILocation(line: 53, column: 1, scope: !1555)
!1590 = !DILocation(line: 54, column: 1, scope: !1555)
!1591 = !DILocation(line: 56, column: 1, scope: !1555)
!1592 = !DILocation(line: 57, column: 1, scope: !1555)
!1593 = !DILocation(line: 59, column: 1, scope: !1555)
!1594 = !DILocation(line: 60, column: 1, scope: !1555)
!1595 = !DILocation(line: 166, column: 1, scope: !1539)
!1596 = !DILocalVariable(scope: !1539, name: "result", file: !3, type: !731)
!1597 = !DILocation(line: 167, column: 1, scope: !1539)
!1598 = !DILocation(line: 168, column: 1, scope: !1539)
!1599 = !DILocalVariable(scope: !1539, name: "k1", file: !3, type: !61)
!1600 = !DILocation(line: 169, column: 1, scope: !1539)
!1601 = !DILocation(line: 72, column: 1, scope: !1539)
!1602 = !DILocation(line: 74, column: 1, scope: !1539)
!1603 = !DILocation(line: 87, column: 1, scope: !1561)
!1604 = !DILocation(line: 76, column: 1, scope: !1559)
!1605 = !DILocation(line: 83, column: 1, scope: !1539)
!1606 = !DILocation(line: 85, column: 1, scope: !1539)
!1607 = !DILocation(line: 86, column: 1, scope: !1539)
!1608 = !DILocation(line: 88, column: 1, scope: !1563)
!1609 = !DILocation(line: 170, column: 1, scope: !1539)
!1610 = !DILocation(line: 172, column: 1, scope: !1539)
!1611 = !DILocation(line: 173, column: 1, scope: !1539)
!1612 = !DILocation(line: 174, column: 1, scope: !1539)
!1613 = !DILocation(line: 176, column: 1, scope: !1539)
!1614 = !DILocation(line: 177, column: 1, scope: !1539)
!1615 = !{ !1669, !1671, !1673, !1675 }
!1616 = distinct !DISubprogram(file: !3, scope: !10, name: "Romberg", line: 180, type: !1535, spFlags: 8, unit: !10, scopeLine: 180)
!1617 = !DILocation(scope: !1616)
!1618 = !DILexicalBlock(file: !3, scope: !1616, line: 180, column: 1)
!1619 = !DILocation(scope: !1618)
!1620 = !DILexicalBlock(file: !3, scope: !1618, line: 1, column: 1)
!1621 = !DILocation(scope: !1620)
!1622 = !DILexicalBlock(file: !3, scope: !1620, line: 52, column: 1)
!1623 = !DILocation(scope: !1622)
!1624 = !DILexicalBlock(file: !3, scope: !1618, line: 1, column: 1)
!1625 = !DILocation(scope: !1624)
!1626 = !DILexicalBlock(file: !3, scope: !1624, line: 76, column: 1)
!1627 = !DILocation(scope: !1626)
!1628 = !DILexicalBlock(file: !3, scope: !1624, line: 87, column: 1)
!1629 = !DILocation(scope: !1628)
!1630 = !DILexicalBlock(file: !3, scope: !1628, line: 88, column: 1)
!1631 = !DILocation(scope: !1630)
!1632 = !DILexicalBlock(file: !3, scope: !1618, line: 1, column: 1)
!1633 = !DILocation(scope: !1632)
!1634 = !DILexicalBlock(file: !3, scope: !1632, line: 1, column: 1)
!1635 = !DILocation(scope: !1634)
!1636 = !DILexicalBlock(file: !3, scope: !1634, line: 1, column: 1)
!1637 = !DILocation(scope: !1636)
!1638 = !DILexicalBlock(file: !3, scope: !1636, line: 1, column: 1)
!1639 = !DILocation(scope: !1638)
!1640 = !DILexicalBlock(file: !3, scope: !1638, line: 1, column: 1)
!1641 = !DILocation(scope: !1640)
!1642 = !DILexicalBlock(file: !3, scope: !1634, line: 1, column: 1)
!1643 = !DILocation(scope: !1642)
!1644 = !DILexicalBlock(file: !3, scope: !1618, line: 1, column: 1)
!1645 = !DILocation(scope: !1644)
!1646 = !DILexicalBlock(file: !3, scope: !1644, line: 52, column: 1)
!1647 = !DILocation(scope: !1646)
!1648 = !DILexicalBlock(file: !3, scope: !1618, line: 1, column: 1)
!1649 = !DILocation(scope: !1648)
!1650 = !DILexicalBlock(file: !3, scope: !1648, line: 76, column: 1)
!1651 = !DILocation(scope: !1650)
!1652 = !DILexicalBlock(file: !3, scope: !1648, line: 87, column: 1)
!1653 = !DILocation(scope: !1652)
!1654 = !DILexicalBlock(file: !3, scope: !1652, line: 88, column: 1)
!1655 = !DILocation(scope: !1654)
!1656 = !DILexicalBlock(file: !3, scope: !1618, line: 1, column: 1)
!1657 = !DILocation(scope: !1656)
!1658 = !DILexicalBlock(file: !3, scope: !1656, line: 1, column: 1)
!1659 = !DILocation(scope: !1658)
!1660 = !DILexicalBlock(file: !3, scope: !1658, line: 1, column: 1)
!1661 = !DILocation(scope: !1660)
!1662 = !DILexicalBlock(file: !3, scope: !1660, line: 1, column: 1)
!1663 = !DILocation(scope: !1662)
!1664 = !DILexicalBlock(file: !3, scope: !1662, line: 1, column: 1)
!1665 = !DILocation(scope: !1664)
!1666 = !DILexicalBlock(file: !3, scope: !1658, line: 1, column: 1)
!1667 = !DILocation(scope: !1666)
!1668 = !DILocalVariable(scope: !1618, name: "func", file: !3, type: !1354)
!1669 = !DILocalVariable(scope: !1616, name: "func", arg: 1, file: !3, type: !1354)
!1670 = !DILocalVariable(scope: !1618, name: "a", file: !3, type: !731)
!1671 = !DILocalVariable(scope: !1616, name: "a", arg: 2, file: !3, type: !731)
!1672 = !DILocalVariable(scope: !1618, name: "b", file: !3, type: !731)
!1673 = !DILocalVariable(scope: !1616, name: "b", arg: 3, file: !3, type: !731)
!1674 = !DILocalVariable(scope: !1618, name: "absacc", file: !3, type: !731)
!1675 = !DILocalVariable(scope: !1616, name: "absacc", arg: 4, file: !3, type: !731)
!1676 = !DILocation(line: 181, column: 1, scope: !1618)
!1677 = !DILocalVariable(scope: !1618, name: "it", file: !3, type: !61)
!1678 = !DILocation(line: 187, column: 1, scope: !1618)
!1679 = !DILocalVariable(scope: !1618, name: "H", file: !3, type: !1578)
!1680 = !DILocation(line: 189, column: 1, scope: !1618)
!1681 = !DILocalVariable(scope: !1618, name: "dresult_pol", file: !3, type: !731)
!1682 = !DILocation(line: 190, column: 1, scope: !1618)
!1683 = !DILocalVariable(scope: !1618, name: "j", file: !3, type: !61)
!1684 = !DILocalVariable(scope: !1618, name: "S", file: !3, type: !1578)
!1685 = !DILocation(line: 191, column: 1, scope: !1618)
!1686 = !DILocation(line: 49, column: 1, scope: !1618)
!1687 = !DILocation(line: 50, column: 1, scope: !1618)
!1688 = !DILocation(line: 52, column: 1, scope: !1646)
!1689 = !DILocation(line: 53, column: 1, scope: !1646)
!1690 = !DILocation(line: 54, column: 1, scope: !1646)
!1691 = !DILocation(line: 56, column: 1, scope: !1646)
!1692 = !DILocation(line: 57, column: 1, scope: !1646)
!1693 = !DILocation(line: 59, column: 1, scope: !1646)
!1694 = !DILocation(line: 60, column: 1, scope: !1646)
!1695 = !DILocation(line: 192, column: 1, scope: !1618)
!1696 = !DILocalVariable(scope: !1618, name: "result", file: !3, type: !731)
!1697 = !DILocation(line: 193, column: 1, scope: !1618)
!1698 = !DILocation(line: 194, column: 1, scope: !1618)
!1699 = !DILocalVariable(scope: !1618, name: "k1", file: !3, type: !61)
!1700 = !DILocation(line: 197, column: 1, scope: !1618)
!1701 = !DILocation(line: 72, column: 1, scope: !1618)
!1702 = !DILocation(line: 74, column: 1, scope: !1618)
!1703 = !DILocation(line: 136, column: 1, scope: !1618)
!1704 = !DILocation(line: 76, column: 1, scope: !1650)
!1705 = !DILocation(line: 83, column: 1, scope: !1618)
!1706 = !DILocalVariable(scope: !1618, name: "result_pol", file: !3, type: !731)
!1707 = !DILocation(line: 85, column: 1, scope: !1618)
!1708 = !DILocation(line: 86, column: 1, scope: !1618)
!1709 = !DILocation(line: 142, column: 1, scope: !1618)
!1710 = !DILocation(line: 87, column: 1, scope: !1652)
!1711 = !DILocation(line: 88, column: 1, scope: !1654)
!1712 = !DILocation(line: 198, column: 1, scope: !1618)
!1713 = !DILocation(line: 116, column: 1, scope: !1618)
!1714 = !DILocation(line: 118, column: 1, scope: !1618)
!1715 = !DILocation(line: 119, column: 1, scope: !1618)
!1716 = !DILocation(line: 120, column: 1, scope: !1618)
!1717 = !DILocation(line: 121, column: 1, scope: !1618)
!1718 = !DILocalVariable(scope: !1618, name: "result_rat", file: !3, type: !731)
!1719 = !DILocation(line: 122, column: 1, scope: !1618)
!1720 = !DILocalVariable(scope: !1618, name: "dresult_rat", file: !3, type: !731)
!1721 = !DILocation(line: 123, column: 1, scope: !1618)
!1722 = !DILocation(line: 125, column: 1, scope: !1618)
!1723 = !DILocation(line: 126, column: 1, scope: !1618)
!1724 = !DILocation(line: 129, column: 1, scope: !1618)
!1725 = !DILocation(line: 130, column: 1, scope: !1618)
!1726 = !DILocation(line: 131, column: 1, scope: !1618)
!1727 = !DILocation(line: 133, column: 1, scope: !1618)
!1728 = !DILocation(line: 135, column: 1, scope: !1618)
!1729 = !DILocation(line: 138, column: 1, scope: !1618)
!1730 = !DILocation(line: 139, column: 1, scope: !1618)
!1731 = !DILocation(line: 140, column: 1, scope: !1618)
!1732 = !DILocation(line: 143, column: 1, scope: !1618)
!1733 = !DILocation(line: 146, column: 1, scope: !1618)
!1734 = !DILocation(line: 147, column: 1, scope: !1618)
!1735 = !DILocation(line: 148, column: 1, scope: !1618)
!1736 = !DILocation(line: 149, column: 1, scope: !1618)
!1737 = !DILocation(line: 151, column: 1, scope: !1618)
!1738 = !DILocation(line: 152, column: 1, scope: !1618)
!1739 = !DILocation(line: 199, column: 1, scope: !1618)
!1740 = !DILocation(line: 200, column: 1, scope: !1618)
!1741 = !DILocation(line: 201, column: 1, scope: !1618)
!1742 = !DILocation(line: 202, column: 1, scope: !1618)
!1743 = !DILocalVariable(scope: !1618, name: "dresult_polrat", file: !3, type: !731)
!1744 = !DILocation(line: 203, column: 1, scope: !1618)
!1745 = !DILocalVariable(scope: !1618, name: "dresult", file: !3, type: !731)
!1746 = !DILocation(line: 204, column: 1, scope: !1618)
!1747 = !DILocation(line: 206, column: 1, scope: !1618)
!1748 = !DILocation(line: 207, column: 1, scope: !1618)
!1749 = !DILocation(line: 210, column: 1, scope: !1618)
!1750 = !DILocation(line: 211, column: 1, scope: !1618)
!1751 = !DILocation(line: 212, column: 1, scope: !1618)
!1752 = !DILocation(line: 214, column: 1, scope: !1618)
!1753 = !DILocation(line: 215, column: 1, scope: !1618)
!1754 = !{ !1756 }
!1755 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "T2DFunction", size: 64, align: 64, elements: !1754)
!1756 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1755, name: "__vptr", size: 64, align: 64, baseType: !126)
!1757 = !DIDerivedType(tag: DW_TAG_reference_type, size: 64, align: 64, baseType: !1755)
!1758 = !{ !731, !1757, !731, !731, !731, !731, !731 }
!1759 = !DISubroutineType(types: !1758)
!1760 = !{ !1854, !1856, !1858, !1860, !1862, !1864 }
!1761 = distinct !DISubprogram(file: !3, scope: !10, name: "Romberg", line: 230, type: !1759, spFlags: 8, unit: !10, scopeLine: 230)
!1762 = !DILocation(scope: !1761)
!1763 = !DILexicalBlock(file: !3, scope: !1761, line: 230, column: 1)
!1764 = !DILocation(scope: !1763)
!1765 = !DILexicalBlock(file: !3, scope: !1763, line: 1, column: 1)
!1766 = !DILocation(scope: !1765)
!1767 = !DILexicalBlock(file: !3, scope: !1765, line: 1, column: 1)
!1768 = !DILocation(scope: !1767)
!1769 = !DILexicalBlock(file: !3, scope: !1767, line: 1, column: 1)
!1770 = !DILocation(scope: !1769)
!1771 = !DILexicalBlock(file: !3, scope: !1763, line: 1, column: 1)
!1772 = !DILocation(scope: !1771)
!1773 = !DILexicalBlock(file: !3, scope: !1771, line: 1, column: 1)
!1774 = !DILocation(scope: !1773)
!1775 = !DILexicalBlock(file: !3, scope: !1773, line: 52, column: 1)
!1776 = !DILocation(scope: !1775)
!1777 = !DILexicalBlock(file: !3, scope: !1771, line: 1, column: 1)
!1778 = !DILocation(scope: !1777)
!1779 = !DILexicalBlock(file: !3, scope: !1777, line: 76, column: 1)
!1780 = !DILocation(scope: !1779)
!1781 = !DILexicalBlock(file: !3, scope: !1777, line: 87, column: 1)
!1782 = !DILocation(scope: !1781)
!1783 = !DILexicalBlock(file: !3, scope: !1781, line: 88, column: 1)
!1784 = !DILocation(scope: !1783)
!1785 = !DILexicalBlock(file: !3, scope: !1771, line: 1, column: 1)
!1786 = !DILocation(scope: !1785)
!1787 = !DILexicalBlock(file: !3, scope: !1785, line: 1, column: 1)
!1788 = !DILocation(scope: !1787)
!1789 = !DILexicalBlock(file: !3, scope: !1787, line: 1, column: 1)
!1790 = !DILocation(scope: !1789)
!1791 = !DILexicalBlock(file: !3, scope: !1789, line: 1, column: 1)
!1792 = !DILocation(scope: !1791)
!1793 = !DILexicalBlock(file: !3, scope: !1791, line: 1, column: 1)
!1794 = !DILocation(scope: !1793)
!1795 = !DILexicalBlock(file: !3, scope: !1787, line: 1, column: 1)
!1796 = !DILocation(scope: !1795)
!1797 = !DILexicalBlock(file: !3, scope: !1763, line: 1, column: 1)
!1798 = !DILocation(scope: !1797)
!1799 = !DILexicalBlock(file: !3, scope: !1797, line: 1, column: 1)
!1800 = !DILocation(scope: !1799)
!1801 = !DILexicalBlock(file: !3, scope: !1799, line: 1, column: 1)
!1802 = !DILocation(scope: !1801)
!1803 = !DILexicalBlock(file: !3, scope: !1763, line: 1, column: 1)
!1804 = !DILocation(scope: !1803)
!1805 = !DILexicalBlock(file: !3, scope: !1803, line: 1, column: 1)
!1806 = !DILocation(scope: !1805)
!1807 = !DILexicalBlock(file: !3, scope: !1805, line: 1, column: 1)
!1808 = !DILocation(scope: !1807)
!1809 = !DILexicalBlock(file: !3, scope: !1763, line: 1, column: 1)
!1810 = !DILocation(scope: !1809)
!1811 = !DILexicalBlock(file: !3, scope: !1809, line: 1, column: 1)
!1812 = !DILocation(scope: !1811)
!1813 = !DILexicalBlock(file: !3, scope: !1811, line: 1, column: 1)
!1814 = !DILocation(scope: !1813)
!1815 = !DILexicalBlock(file: !3, scope: !1763, line: 1, column: 1)
!1816 = !DILocation(scope: !1815)
!1817 = !DILexicalBlock(file: !3, scope: !1815, line: 1, column: 1)
!1818 = !DILocation(scope: !1817)
!1819 = !DILexicalBlock(file: !3, scope: !1817, line: 52, column: 1)
!1820 = !DILocation(scope: !1819)
!1821 = !DILexicalBlock(file: !3, scope: !1815, line: 1, column: 1)
!1822 = !DILocation(scope: !1821)
!1823 = !DILexicalBlock(file: !3, scope: !1821, line: 76, column: 1)
!1824 = !DILocation(scope: !1823)
!1825 = !DILexicalBlock(file: !3, scope: !1821, line: 87, column: 1)
!1826 = !DILocation(scope: !1825)
!1827 = !DILexicalBlock(file: !3, scope: !1825, line: 88, column: 1)
!1828 = !DILocation(scope: !1827)
!1829 = !DILexicalBlock(file: !3, scope: !1815, line: 1, column: 1)
!1830 = !DILocation(scope: !1829)
!1831 = !DILexicalBlock(file: !3, scope: !1829, line: 1, column: 1)
!1832 = !DILocation(scope: !1831)
!1833 = !DILexicalBlock(file: !3, scope: !1831, line: 1, column: 1)
!1834 = !DILocation(scope: !1833)
!1835 = !DILexicalBlock(file: !3, scope: !1833, line: 1, column: 1)
!1836 = !DILocation(scope: !1835)
!1837 = !DILexicalBlock(file: !3, scope: !1835, line: 1, column: 1)
!1838 = !DILocation(scope: !1837)
!1839 = !DILexicalBlock(file: !3, scope: !1831, line: 1, column: 1)
!1840 = !DILocation(scope: !1839)
!1841 = !DILexicalBlock(file: !3, scope: !1763, line: 1, column: 1)
!1842 = !DILocation(scope: !1841)
!1843 = !DILexicalBlock(file: !3, scope: !1841, line: 1, column: 1)
!1844 = !DILocation(scope: !1843)
!1845 = !DILexicalBlock(file: !3, scope: !1843, line: 1, column: 1)
!1846 = !DILocation(scope: !1845)
!1847 = !DILexicalBlock(file: !3, scope: !1763, line: 1, column: 1)
!1848 = !DILocation(scope: !1847)
!1849 = !DILexicalBlock(file: !3, scope: !1847, line: 1, column: 1)
!1850 = !DILocation(scope: !1849)
!1851 = !DILexicalBlock(file: !3, scope: !1849, line: 1, column: 1)
!1852 = !DILocation(scope: !1851)
!1853 = !DILocalVariable(scope: !1763, name: "func", file: !3, type: !1757)
!1854 = !DILocalVariable(scope: !1761, name: "func", arg: 1, file: !3, type: !1757)
!1855 = !DILocalVariable(scope: !1763, name: "a", file: !3, type: !731)
!1856 = !DILocalVariable(scope: !1761, name: "a", arg: 2, file: !3, type: !731)
!1857 = !DILocalVariable(scope: !1763, name: "b", file: !3, type: !731)
!1858 = !DILocalVariable(scope: !1761, name: "b", arg: 3, file: !3, type: !731)
!1859 = !DILocalVariable(scope: !1763, name: "c", file: !3, type: !731)
!1860 = !DILocalVariable(scope: !1761, name: "c", arg: 4, file: !3, type: !731)
!1861 = !DILocalVariable(scope: !1763, name: "d", file: !3, type: !731)
!1862 = !DILocalVariable(scope: !1761, name: "d", arg: 5, file: !3, type: !731)
!1863 = !DILocalVariable(scope: !1763, name: "absacc", file: !3, type: !731)
!1864 = !DILocalVariable(scope: !1761, name: "absacc", arg: 6, file: !3, type: !731)
!1865 = !DILocation(line: 232, column: 1, scope: !1763)
!1866 = !DILocation(line: 225, column: 1, scope: !1763)
!1867 = !DISubrange(count: 5)
!1868 = !{ !1867 }
!1869 = !DICompositeType(tag: DW_TAG_array_type, size: 320, align: 64, baseType: !126, elements: !1868)
!1870 = distinct !DIGlobalVariable(scope: !10, name: "_ZTV11T1DFunction", file: !3, type: !1869, isDefinition: true)
!1871 = !DIGlobalVariableExpression(var: !1870, expr: !1366)
!1872 = distinct !DIGlobalVariable(scope: !10, name: "_ZTV9Tinty_f2D", file: !3, type: !1869, isDefinition: true)
!1873 = !DIGlobalVariableExpression(var: !1872, expr: !1366)
!1874 = !DILocation(line: 184, column: 1, scope: !1763)
!1875 = !DILocation(line: 188, column: 1, scope: !1763)
!1876 = !DILocation(line: 189, column: 1, scope: !1763)
!1877 = !DILocation(line: 190, column: 1, scope: !1763)
!1878 = !DILocation(line: 191, column: 1, scope: !1763)
!1879 = !DILocation(line: 192, column: 1, scope: !1763)
!1880 = !DILocation(line: 194, column: 1, scope: !1763)
!1881 = !DILocation(line: 197, column: 1, scope: !1763)
!1882 = !DILocation(line: 198, column: 1, scope: !1763)
!1883 = !DILocation(line: 199, column: 1, scope: !1763)
!1884 = !DILocation(line: 200, column: 1, scope: !1763)
!1885 = !DILocation(line: 201, column: 1, scope: !1763)
!1886 = !DILocation(line: 202, column: 1, scope: !1763)
!1887 = !DILocation(line: 204, column: 1, scope: !1763)
!1888 = !DILocation(line: 206, column: 1, scope: !1763)
!1889 = !DILocation(line: 210, column: 1, scope: !1763)
!1890 = !DILocation(line: 211, column: 1, scope: !1763)
!1891 = !DILocation(line: 212, column: 1, scope: !1763)
!1892 = !DILocation(line: 215, column: 1, scope: !1763)
!1893 = !DILocation(line: 233, column: 1, scope: !1763)
!1894 = !{ !1896 }
!1895 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "T3DFunction", size: 64, align: 64, elements: !1894)
!1896 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !1895, name: "__vptr", size: 64, align: 64, baseType: !126)
!1897 = !DIDerivedType(tag: DW_TAG_reference_type, size: 64, align: 64, baseType: !1895)
!1898 = !{ !731, !1897, !731, !731, !731, !731, !731, !731, !731 }
!1899 = !DISubroutineType(types: !1898)
!1900 = !{ !1994, !1996, !1998, !2000, !2002, !2004, !2006, !2008 }
!1901 = distinct !DISubprogram(file: !3, scope: !10, name: "Romberg", line: 248, type: !1899, spFlags: 8, unit: !10, scopeLine: 248)
!1902 = !DILocation(scope: !1901)
!1903 = !DILexicalBlock(file: !3, scope: !1901, line: 248, column: 1)
!1904 = !DILocation(scope: !1903)
!1905 = !DILexicalBlock(file: !3, scope: !1903, line: 1, column: 1)
!1906 = !DILocation(scope: !1905)
!1907 = !DILexicalBlock(file: !3, scope: !1905, line: 1, column: 1)
!1908 = !DILocation(scope: !1907)
!1909 = !DILexicalBlock(file: !3, scope: !1907, line: 1, column: 1)
!1910 = !DILocation(scope: !1909)
!1911 = !DILexicalBlock(file: !3, scope: !1903, line: 1, column: 1)
!1912 = !DILocation(scope: !1911)
!1913 = !DILexicalBlock(file: !3, scope: !1911, line: 1, column: 1)
!1914 = !DILocation(scope: !1913)
!1915 = !DILexicalBlock(file: !3, scope: !1913, line: 52, column: 1)
!1916 = !DILocation(scope: !1915)
!1917 = !DILexicalBlock(file: !3, scope: !1911, line: 1, column: 1)
!1918 = !DILocation(scope: !1917)
!1919 = !DILexicalBlock(file: !3, scope: !1917, line: 76, column: 1)
!1920 = !DILocation(scope: !1919)
!1921 = !DILexicalBlock(file: !3, scope: !1917, line: 87, column: 1)
!1922 = !DILocation(scope: !1921)
!1923 = !DILexicalBlock(file: !3, scope: !1921, line: 88, column: 1)
!1924 = !DILocation(scope: !1923)
!1925 = !DILexicalBlock(file: !3, scope: !1911, line: 1, column: 1)
!1926 = !DILocation(scope: !1925)
!1927 = !DILexicalBlock(file: !3, scope: !1925, line: 1, column: 1)
!1928 = !DILocation(scope: !1927)
!1929 = !DILexicalBlock(file: !3, scope: !1927, line: 1, column: 1)
!1930 = !DILocation(scope: !1929)
!1931 = !DILexicalBlock(file: !3, scope: !1929, line: 1, column: 1)
!1932 = !DILocation(scope: !1931)
!1933 = !DILexicalBlock(file: !3, scope: !1931, line: 1, column: 1)
!1934 = !DILocation(scope: !1933)
!1935 = !DILexicalBlock(file: !3, scope: !1927, line: 1, column: 1)
!1936 = !DILocation(scope: !1935)
!1937 = !DILexicalBlock(file: !3, scope: !1903, line: 1, column: 1)
!1938 = !DILocation(scope: !1937)
!1939 = !DILexicalBlock(file: !3, scope: !1937, line: 1, column: 1)
!1940 = !DILocation(scope: !1939)
!1941 = !DILexicalBlock(file: !3, scope: !1939, line: 1, column: 1)
!1942 = !DILocation(scope: !1941)
!1943 = !DILexicalBlock(file: !3, scope: !1903, line: 1, column: 1)
!1944 = !DILocation(scope: !1943)
!1945 = !DILexicalBlock(file: !3, scope: !1943, line: 1, column: 1)
!1946 = !DILocation(scope: !1945)
!1947 = !DILexicalBlock(file: !3, scope: !1945, line: 1, column: 1)
!1948 = !DILocation(scope: !1947)
!1949 = !DILexicalBlock(file: !3, scope: !1903, line: 1, column: 1)
!1950 = !DILocation(scope: !1949)
!1951 = !DILexicalBlock(file: !3, scope: !1949, line: 1, column: 1)
!1952 = !DILocation(scope: !1951)
!1953 = !DILexicalBlock(file: !3, scope: !1951, line: 1, column: 1)
!1954 = !DILocation(scope: !1953)
!1955 = !DILexicalBlock(file: !3, scope: !1903, line: 1, column: 1)
!1956 = !DILocation(scope: !1955)
!1957 = !DILexicalBlock(file: !3, scope: !1955, line: 1, column: 1)
!1958 = !DILocation(scope: !1957)
!1959 = !DILexicalBlock(file: !3, scope: !1957, line: 52, column: 1)
!1960 = !DILocation(scope: !1959)
!1961 = !DILexicalBlock(file: !3, scope: !1955, line: 1, column: 1)
!1962 = !DILocation(scope: !1961)
!1963 = !DILexicalBlock(file: !3, scope: !1961, line: 76, column: 1)
!1964 = !DILocation(scope: !1963)
!1965 = !DILexicalBlock(file: !3, scope: !1961, line: 87, column: 1)
!1966 = !DILocation(scope: !1965)
!1967 = !DILexicalBlock(file: !3, scope: !1965, line: 88, column: 1)
!1968 = !DILocation(scope: !1967)
!1969 = !DILexicalBlock(file: !3, scope: !1955, line: 1, column: 1)
!1970 = !DILocation(scope: !1969)
!1971 = !DILexicalBlock(file: !3, scope: !1969, line: 1, column: 1)
!1972 = !DILocation(scope: !1971)
!1973 = !DILexicalBlock(file: !3, scope: !1971, line: 1, column: 1)
!1974 = !DILocation(scope: !1973)
!1975 = !DILexicalBlock(file: !3, scope: !1973, line: 1, column: 1)
!1976 = !DILocation(scope: !1975)
!1977 = !DILexicalBlock(file: !3, scope: !1975, line: 1, column: 1)
!1978 = !DILocation(scope: !1977)
!1979 = !DILexicalBlock(file: !3, scope: !1971, line: 1, column: 1)
!1980 = !DILocation(scope: !1979)
!1981 = !DILexicalBlock(file: !3, scope: !1903, line: 1, column: 1)
!1982 = !DILocation(scope: !1981)
!1983 = !DILexicalBlock(file: !3, scope: !1981, line: 1, column: 1)
!1984 = !DILocation(scope: !1983)
!1985 = !DILexicalBlock(file: !3, scope: !1983, line: 1, column: 1)
!1986 = !DILocation(scope: !1985)
!1987 = !DILexicalBlock(file: !3, scope: !1903, line: 1, column: 1)
!1988 = !DILocation(scope: !1987)
!1989 = !DILexicalBlock(file: !3, scope: !1987, line: 1, column: 1)
!1990 = !DILocation(scope: !1989)
!1991 = !DILexicalBlock(file: !3, scope: !1989, line: 1, column: 1)
!1992 = !DILocation(scope: !1991)
!1993 = !DILocalVariable(scope: !1903, name: "func", file: !3, type: !1897)
!1994 = !DILocalVariable(scope: !1901, name: "func", arg: 1, file: !3, type: !1897)
!1995 = !DILocalVariable(scope: !1903, name: "a", file: !3, type: !731)
!1996 = !DILocalVariable(scope: !1901, name: "a", arg: 2, file: !3, type: !731)
!1997 = !DILocalVariable(scope: !1903, name: "b", file: !3, type: !731)
!1998 = !DILocalVariable(scope: !1901, name: "b", arg: 3, file: !3, type: !731)
!1999 = !DILocalVariable(scope: !1903, name: "c", file: !3, type: !731)
!2000 = !DILocalVariable(scope: !1901, name: "c", arg: 4, file: !3, type: !731)
!2001 = !DILocalVariable(scope: !1903, name: "d", file: !3, type: !731)
!2002 = !DILocalVariable(scope: !1901, name: "d", arg: 5, file: !3, type: !731)
!2003 = !DILocalVariable(scope: !1903, name: "e", file: !3, type: !731)
!2004 = !DILocalVariable(scope: !1901, name: "e", arg: 6, file: !3, type: !731)
!2005 = !DILocalVariable(scope: !1903, name: "f", file: !3, type: !731)
!2006 = !DILocalVariable(scope: !1901, name: "f", arg: 7, file: !3, type: !731)
!2007 = !DILocalVariable(scope: !1903, name: "absacc", file: !3, type: !731)
!2008 = !DILocalVariable(scope: !1901, name: "absacc", arg: 8, file: !3, type: !731)
!2009 = !DILocation(line: 249, column: 1, scope: !1903)
!2010 = !DILocation(line: 243, column: 1, scope: !1903)
!2011 = distinct !DIGlobalVariable(scope: !10, name: "_ZTV10Tintxy_f3D", file: !3, type: !1869, isDefinition: true)
!2012 = !DIGlobalVariableExpression(var: !2011, expr: !1366)
!2013 = !DILocation(line: 184, column: 1, scope: !1903)
!2014 = !DILocation(line: 188, column: 1, scope: !1903)
!2015 = !DILocation(line: 189, column: 1, scope: !1903)
!2016 = !DILocation(line: 190, column: 1, scope: !1903)
!2017 = !DILocation(line: 191, column: 1, scope: !1903)
!2018 = !DILocation(line: 192, column: 1, scope: !1903)
!2019 = !DILocation(line: 194, column: 1, scope: !1903)
!2020 = !DILocation(line: 197, column: 1, scope: !1903)
!2021 = !DILocation(line: 198, column: 1, scope: !1903)
!2022 = !DILocation(line: 199, column: 1, scope: !1903)
!2023 = !DILocation(line: 200, column: 1, scope: !1903)
!2024 = !DILocation(line: 201, column: 1, scope: !1903)
!2025 = !DILocation(line: 202, column: 1, scope: !1903)
!2026 = !DILocation(line: 204, column: 1, scope: !1903)
!2027 = !DILocation(line: 206, column: 1, scope: !1903)
!2028 = !DILocation(line: 210, column: 1, scope: !1903)
!2029 = !DILocation(line: 211, column: 1, scope: !1903)
!2030 = !DILocation(line: 212, column: 1, scope: !1903)
!2031 = !DILocation(line: 215, column: 1, scope: !1903)
!2032 = !DILocation(line: 250, column: 1, scope: !1903)
!2033 = !DIFile(filename: "backgroundfield/functions.hpp", directory: "/home/talgat/vlasiator")
; !2034 = !DIFile(tag: DW_TAG_file_type, pair: !2033)
!2034 = !{ i32 41, !2033 }
!2035 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1352)
!2036 = !{ null, !2035 }
!2037 = !DISubroutineType(types: !2036)
!2038 = !{ !2044 }
!2039 = distinct !DISubprogram(file: !2033, scope: !1352, name: "~T1DFunction", line: 29, type: !2037, spFlags: 8, unit: !10, scopeLine: 29)
!2040 = !DILocation(scope: !2039)
!2041 = !DILexicalBlock(file: !2033, scope: !2039, line: 29, column: 1)
!2042 = !DILocation(scope: !2041)
!2043 = !DILocalVariable(scope: !2041, file: !2033, type: !2035, flags: 64)
!2044 = !DILocalVariable(scope: !2039, arg: 1, file: !2033, type: !2035, flags: 64)
!2045 = !DILocation(line: 29, column: 1, scope: !2041)
!2046 = !{  }
!2047 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN11T1DFunctionD0Ev", type: !2037, spFlags: 8, unit: !10)
!2048 = !DILocation(scope: !2047)
!2049 = !DILexicalBlock(file: !3, scope: !2047, line: 1, column: 1)
!2050 = !DILocation(scope: !2049)
!2051 = !DILexicalBlock(file: !3, scope: !2049, line: 1, column: 1)
!2052 = !DILocation(scope: !2051)
!2053 = !DILexicalBlock(file: !3, scope: !2049, line: 1, column: 1)
!2054 = !DILocation(scope: !2053)
!2055 = !DILocation(line: 29, column: 1, scope: !2049)
!2056 = !{  }
!2057 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN11T1DFunctionD2Ev", type: !2037, spFlags: 8, unit: !10)
!2058 = !DILocation(scope: !2057)
!2059 = !DILexicalBlock(file: !3, scope: !2057, line: 1, column: 1)
!2060 = !DILocation(scope: !2059)
!2061 = !DILexicalBlock(file: !3, scope: !2059, line: 1, column: 1)
!2062 = !DILocation(scope: !2061)
!2063 = !DILexicalBlock(file: !3, scope: !2059, line: 1, column: 1)
!2064 = !DILocation(scope: !2063)
!2065 = !DILocation(line: 29, column: 1, scope: !2059)
!2066 = !{ !2072 }
!2067 = distinct !DISubprogram(file: !2033, scope: !10, name: "T1DFunction", type: !2037, spFlags: 8, unit: !10)
!2068 = !DILocation(scope: !2067)
!2069 = !DILexicalBlock(file: !2033, scope: !2067, line: 1, column: 1)
!2070 = !DILocation(scope: !2069)
!2071 = !DILocalVariable(scope: !2069, file: !2033, type: !2035, flags: 64)
!2072 = !DILocalVariable(scope: !2067, arg: 1, file: !2033, type: !2035, flags: 64)
!2073 = !DILocation(line: 29, column: 1, scope: !2069)
!2074 = !{  }
!2075 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN11T1DFunctionC2Ev", type: !2037, spFlags: 8, unit: !10)
!2076 = !DILocation(scope: !2075)
!2077 = !DILexicalBlock(file: !3, scope: !2075, line: 1, column: 1)
!2078 = !DILocation(scope: !2077)
!2079 = !DILexicalBlock(file: !3, scope: !2077, line: 1, column: 1)
!2080 = !DILocation(scope: !2079)
!2081 = !DILexicalBlock(file: !3, scope: !2077, line: 1, column: 1)
!2082 = !DILocation(scope: !2081)
!2083 = !DILocation(line: 29, column: 1, scope: !2077)
!2084 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !1755)
!2085 = !{ null, !2084 }
!2086 = !DISubroutineType(types: !2085)
!2087 = !{ !2093 }
!2088 = distinct !DISubprogram(file: !2033, scope: !1755, name: "~T2DFunction", line: 30, type: !2086, spFlags: 8, unit: !10, scopeLine: 30)
!2089 = !DILocation(scope: !2088)
!2090 = !DILexicalBlock(file: !2033, scope: !2088, line: 30, column: 1)
!2091 = !DILocation(scope: !2090)
!2092 = !DILocalVariable(scope: !2090, file: !2033, type: !2084, flags: 64)
!2093 = !DILocalVariable(scope: !2088, arg: 1, file: !2033, type: !2084, flags: 64)
!2094 = !DILocation(line: 30, column: 1, scope: !2090)
!2095 = distinct !DIGlobalVariable(scope: !10, name: "_ZTV11T2DFunction", file: !3, type: !1869, isDefinition: true)
!2096 = !DIGlobalVariableExpression(var: !2095, expr: !1366)
!2097 = !{  }
!2098 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN11T2DFunctionD0Ev", type: !2086, spFlags: 8, unit: !10)
!2099 = !DILocation(scope: !2098)
!2100 = !DILexicalBlock(file: !3, scope: !2098, line: 1, column: 1)
!2101 = !DILocation(scope: !2100)
!2102 = !DILexicalBlock(file: !3, scope: !2100, line: 1, column: 1)
!2103 = !DILocation(scope: !2102)
!2104 = !DILexicalBlock(file: !3, scope: !2100, line: 1, column: 1)
!2105 = !DILocation(scope: !2104)
!2106 = !DILocation(line: 30, column: 1, scope: !2100)
!2107 = !{  }
!2108 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN11T2DFunctionD2Ev", type: !2086, spFlags: 8, unit: !10)
!2109 = !DILocation(scope: !2108)
!2110 = !DILexicalBlock(file: !3, scope: !2108, line: 1, column: 1)
!2111 = !DILocation(scope: !2110)
!2112 = !DILexicalBlock(file: !3, scope: !2110, line: 1, column: 1)
!2113 = !DILocation(scope: !2112)
!2114 = !DILexicalBlock(file: !3, scope: !2110, line: 1, column: 1)
!2115 = !DILocation(scope: !2114)
!2116 = !DILocation(line: 30, column: 1, scope: !2110)
!2117 = !{ !2123 }
!2118 = distinct !DISubprogram(file: !2033, scope: !10, name: "T2DFunction", type: !2086, spFlags: 8, unit: !10)
!2119 = !DILocation(scope: !2118)
!2120 = !DILexicalBlock(file: !2033, scope: !2118, line: 1, column: 1)
!2121 = !DILocation(scope: !2120)
!2122 = !DILocalVariable(scope: !2120, file: !2033, type: !2084, flags: 64)
!2123 = !DILocalVariable(scope: !2118, arg: 1, file: !2033, type: !2084, flags: 64)
!2124 = !DILocation(line: 30, column: 1, scope: !2120)
!2125 = !{  }
!2126 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN11T2DFunctionC2Ev", type: !2086, spFlags: 8, unit: !10)
!2127 = !DILocation(scope: !2126)
!2128 = !DILexicalBlock(file: !3, scope: !2126, line: 1, column: 1)
!2129 = !DILocation(scope: !2128)
!2130 = !DILexicalBlock(file: !3, scope: !2128, line: 1, column: 1)
!2131 = !DILocation(scope: !2130)
!2132 = !DILexicalBlock(file: !3, scope: !2128, line: 1, column: 1)
!2133 = !DILocation(scope: !2132)
!2134 = !DILocation(line: 30, column: 1, scope: !2128)
!2135 = !{ null, !1188, !1757, !731 }
!2136 = !DISubroutineType(types: !2135)
!2137 = !{ !2151, !2153, !2155 }
!2138 = distinct !DISubprogram(file: !2033, scope: !1183, name: "T2D_fix1", line: 40, type: !2136, spFlags: 8, unit: !10, scopeLine: 40)
!2139 = !DILocation(scope: !2138)
!2140 = !DILexicalBlock(file: !2033, scope: !2138, line: 40, column: 1)
!2141 = !DILocation(scope: !2140)
!2142 = !DILexicalBlock(file: !2033, scope: !2140, line: 1, column: 1)
!2143 = !DILocation(scope: !2142)
!2144 = !DILexicalBlock(file: !2033, scope: !2142, line: 1, column: 1)
!2145 = !DILocation(scope: !2144)
!2146 = !DILexicalBlock(file: !2033, scope: !2140, line: 1, column: 1)
!2147 = !DILocation(scope: !2146)
!2148 = !DILexicalBlock(file: !2033, scope: !2146, line: 1, column: 1)
!2149 = !DILocation(scope: !2148)
!2150 = !DILocalVariable(scope: !2140, file: !2033, type: !1188, flags: 64)
!2151 = !DILocalVariable(scope: !2138, arg: 1, file: !2033, type: !1188, flags: 64)
!2152 = !DILocalVariable(scope: !2140, name: "f1", file: !2033, type: !1757)
!2153 = !DILocalVariable(scope: !2138, name: "f1", arg: 2, file: !2033, type: !1757)
!2154 = !DILocalVariable(scope: !2140, name: "x1", file: !2033, type: !731)
!2155 = !DILocalVariable(scope: !2138, name: "x1", arg: 3, file: !2033, type: !731)
!2156 = !DILocation(line: 40, column: 1, scope: !2140)
!2157 = distinct !DIGlobalVariable(scope: !10, name: "_ZTV8T2D_fix1", file: !3, type: !1869, isDefinition: true)
!2158 = !DIGlobalVariableExpression(var: !2157, expr: !1366)
!2159 = !{  }
!2160 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN8T2D_fix1C2ERK11T2DFunctiond", type: !2136, spFlags: 8, unit: !10)
!2161 = !DILocation(scope: !2160)
!2162 = !DILexicalBlock(file: !3, scope: !2160, line: 1, column: 1)
!2163 = !DILocation(scope: !2162)
!2164 = !DILexicalBlock(file: !3, scope: !2162, line: 1, column: 1)
!2165 = !DILocation(scope: !2164)
!2166 = !DILexicalBlock(file: !3, scope: !2164, line: 1, column: 1)
!2167 = !DILocation(scope: !2166)
!2168 = !DILexicalBlock(file: !3, scope: !2166, line: 1, column: 1)
!2169 = !DILocation(scope: !2168)
!2170 = !DILexicalBlock(file: !3, scope: !2162, line: 1, column: 1)
!2171 = !DILocation(scope: !2170)
!2172 = !DILexicalBlock(file: !3, scope: !2170, line: 1, column: 1)
!2173 = !DILocation(scope: !2172)
!2174 = !DILexicalBlock(file: !3, scope: !2172, line: 1, column: 1)
!2175 = !DILocation(scope: !2174)
!2176 = !DILocation(line: 40, column: 1, scope: !2162)
!2177 = !{ !2183, !2185 }
!2178 = distinct !DISubprogram(file: !2033, scope: !10, name: "call", line: 41, type: !1190, spFlags: 8, unit: !10, scopeLine: 41)
!2179 = !DILocation(scope: !2178)
!2180 = !DILexicalBlock(file: !2033, scope: !2178, line: 41, column: 1)
!2181 = !DILocation(scope: !2180)
!2182 = !DILocalVariable(scope: !2180, file: !2033, type: !1188, flags: 64)
!2183 = !DILocalVariable(scope: !2178, arg: 1, file: !2033, type: !1188, flags: 64)
!2184 = !DILocalVariable(scope: !2180, name: "y", file: !2033, type: !731)
!2185 = !DILocalVariable(scope: !2178, name: "y", arg: 2, file: !2033, type: !731)
!2186 = !DILocation(line: 41, column: 1, scope: !2180)
!2187 = !{ null, !1188 }
!2188 = !DISubroutineType(types: !2187)
!2189 = !{ !2203 }
!2190 = distinct !DISubprogram(file: !2033, scope: !1183, name: "~T2D_fix1", line: 42, type: !2188, spFlags: 8, unit: !10, scopeLine: 42)
!2191 = !DILocation(scope: !2190)
!2192 = !DILexicalBlock(file: !2033, scope: !2190, line: 42, column: 1)
!2193 = !DILocation(scope: !2192)
!2194 = !DILexicalBlock(file: !2033, scope: !2192, line: 1, column: 1)
!2195 = !DILocation(scope: !2194)
!2196 = !DILexicalBlock(file: !2033, scope: !2194, line: 1, column: 1)
!2197 = !DILocation(scope: !2196)
!2198 = !DILexicalBlock(file: !2033, scope: !2192, line: 1, column: 1)
!2199 = !DILocation(scope: !2198)
!2200 = !DILexicalBlock(file: !2033, scope: !2198, line: 1, column: 1)
!2201 = !DILocation(scope: !2200)
!2202 = !DILocalVariable(scope: !2192, file: !2033, type: !1188, flags: 64)
!2203 = !DILocalVariable(scope: !2190, arg: 1, file: !2033, type: !1188, flags: 64)
!2204 = !DILocation(line: 42, column: 1, scope: !2192)
!2205 = !{  }
!2206 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN8T2D_fix1D0Ev", type: !2188, spFlags: 8, unit: !10)
!2207 = !DILocation(scope: !2206)
!2208 = !DILexicalBlock(file: !3, scope: !2206, line: 1, column: 1)
!2209 = !DILocation(scope: !2208)
!2210 = !DILexicalBlock(file: !3, scope: !2208, line: 1, column: 1)
!2211 = !DILocation(scope: !2210)
!2212 = !DILexicalBlock(file: !3, scope: !2210, line: 1, column: 1)
!2213 = !DILocation(scope: !2212)
!2214 = !DILexicalBlock(file: !3, scope: !2212, line: 1, column: 1)
!2215 = !DILocation(scope: !2214)
!2216 = !DILexicalBlock(file: !3, scope: !2208, line: 1, column: 1)
!2217 = !DILocation(scope: !2216)
!2218 = !DILexicalBlock(file: !3, scope: !2216, line: 1, column: 1)
!2219 = !DILocation(scope: !2218)
!2220 = !DILexicalBlock(file: !3, scope: !2218, line: 1, column: 1)
!2221 = !DILocation(scope: !2220)
!2222 = !DILocation(line: 42, column: 1, scope: !2208)
!2223 = !{  }
!2224 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN8T2D_fix1D2Ev", type: !2188, spFlags: 8, unit: !10)
!2225 = !DILocation(scope: !2224)
!2226 = !DILexicalBlock(file: !3, scope: !2224, line: 1, column: 1)
!2227 = !DILocation(scope: !2226)
!2228 = !DILexicalBlock(file: !3, scope: !2226, line: 1, column: 1)
!2229 = !DILocation(scope: !2228)
!2230 = !DILexicalBlock(file: !3, scope: !2228, line: 1, column: 1)
!2231 = !DILocation(scope: !2230)
!2232 = !DILexicalBlock(file: !3, scope: !2230, line: 1, column: 1)
!2233 = !DILocation(scope: !2232)
!2234 = !DILexicalBlock(file: !3, scope: !2226, line: 1, column: 1)
!2235 = !DILocation(scope: !2234)
!2236 = !DILexicalBlock(file: !3, scope: !2234, line: 1, column: 1)
!2237 = !DILocation(scope: !2236)
!2238 = !DILexicalBlock(file: !3, scope: !2236, line: 1, column: 1)
!2239 = !DILocation(scope: !2238)
!2240 = !DILocation(line: 42, column: 1, scope: !2226)
!2241 = !{ null, !1213, !1897, !731 }
!2242 = !DISubroutineType(types: !2241)
!2243 = !{ !2257, !2259, !2261 }
!2244 = distinct !DISubprogram(file: !2033, scope: !1209, name: "T3D_fix3", line: 82, type: !2242, spFlags: 8, unit: !10, scopeLine: 82)
!2245 = !DILocation(scope: !2244)
!2246 = !DILexicalBlock(file: !2033, scope: !2244, line: 82, column: 1)
!2247 = !DILocation(scope: !2246)
!2248 = !DILexicalBlock(file: !2033, scope: !2246, line: 1, column: 1)
!2249 = !DILocation(scope: !2248)
!2250 = !DILexicalBlock(file: !2033, scope: !2248, line: 1, column: 1)
!2251 = !DILocation(scope: !2250)
!2252 = !DILexicalBlock(file: !2033, scope: !2246, line: 1, column: 1)
!2253 = !DILocation(scope: !2252)
!2254 = !DILexicalBlock(file: !2033, scope: !2252, line: 1, column: 1)
!2255 = !DILocation(scope: !2254)
!2256 = !DILocalVariable(scope: !2246, file: !2033, type: !1213, flags: 64)
!2257 = !DILocalVariable(scope: !2244, arg: 1, file: !2033, type: !1213, flags: 64)
!2258 = !DILocalVariable(scope: !2246, name: "f1", file: !2033, type: !1897)
!2259 = !DILocalVariable(scope: !2244, name: "f1", arg: 2, file: !2033, type: !1897)
!2260 = !DILocalVariable(scope: !2246, name: "z1", file: !2033, type: !731)
!2261 = !DILocalVariable(scope: !2244, name: "z1", arg: 3, file: !2033, type: !731)
!2262 = !DILocation(line: 82, column: 1, scope: !2246)
!2263 = distinct !DIGlobalVariable(scope: !10, name: "_ZTV8T3D_fix3", file: !3, type: !1869, isDefinition: true)
!2264 = !DIGlobalVariableExpression(var: !2263, expr: !1366)
!2265 = !{  }
!2266 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN8T3D_fix3C2ERK11T3DFunctiond", type: !2242, spFlags: 8, unit: !10)
!2267 = !DILocation(scope: !2266)
!2268 = !DILexicalBlock(file: !3, scope: !2266, line: 1, column: 1)
!2269 = !DILocation(scope: !2268)
!2270 = !DILexicalBlock(file: !3, scope: !2268, line: 1, column: 1)
!2271 = !DILocation(scope: !2270)
!2272 = !DILexicalBlock(file: !3, scope: !2270, line: 1, column: 1)
!2273 = !DILocation(scope: !2272)
!2274 = !DILexicalBlock(file: !3, scope: !2272, line: 1, column: 1)
!2275 = !DILocation(scope: !2274)
!2276 = !DILexicalBlock(file: !3, scope: !2268, line: 1, column: 1)
!2277 = !DILocation(scope: !2276)
!2278 = !DILexicalBlock(file: !3, scope: !2276, line: 1, column: 1)
!2279 = !DILocation(scope: !2278)
!2280 = !DILexicalBlock(file: !3, scope: !2278, line: 1, column: 1)
!2281 = !DILocation(scope: !2280)
!2282 = !DILocation(line: 82, column: 1, scope: !2268)
!2283 = !{ !2289, !2291, !2293 }
!2284 = distinct !DISubprogram(file: !2033, scope: !10, name: "call", line: 83, type: !1215, spFlags: 8, unit: !10, scopeLine: 83)
!2285 = !DILocation(scope: !2284)
!2286 = !DILexicalBlock(file: !2033, scope: !2284, line: 83, column: 1)
!2287 = !DILocation(scope: !2286)
!2288 = !DILocalVariable(scope: !2286, file: !2033, type: !1213, flags: 64)
!2289 = !DILocalVariable(scope: !2284, arg: 1, file: !2033, type: !1213, flags: 64)
!2290 = !DILocalVariable(scope: !2286, name: "x", file: !2033, type: !731)
!2291 = !DILocalVariable(scope: !2284, name: "x", arg: 2, file: !2033, type: !731)
!2292 = !DILocalVariable(scope: !2286, name: "y", file: !2033, type: !731)
!2293 = !DILocalVariable(scope: !2284, name: "y", arg: 3, file: !2033, type: !731)
!2294 = !DILocation(line: 83, column: 1, scope: !2286)
!2295 = !{ null, !1213 }
!2296 = !DISubroutineType(types: !2295)
!2297 = !{ !2311 }
!2298 = distinct !DISubprogram(file: !2033, scope: !1209, name: "~T3D_fix3", line: 84, type: !2296, spFlags: 8, unit: !10, scopeLine: 84)
!2299 = !DILocation(scope: !2298)
!2300 = !DILexicalBlock(file: !2033, scope: !2298, line: 84, column: 1)
!2301 = !DILocation(scope: !2300)
!2302 = !DILexicalBlock(file: !2033, scope: !2300, line: 1, column: 1)
!2303 = !DILocation(scope: !2302)
!2304 = !DILexicalBlock(file: !2033, scope: !2302, line: 1, column: 1)
!2305 = !DILocation(scope: !2304)
!2306 = !DILexicalBlock(file: !2033, scope: !2300, line: 1, column: 1)
!2307 = !DILocation(scope: !2306)
!2308 = !DILexicalBlock(file: !2033, scope: !2306, line: 1, column: 1)
!2309 = !DILocation(scope: !2308)
!2310 = !DILocalVariable(scope: !2300, file: !2033, type: !1213, flags: 64)
!2311 = !DILocalVariable(scope: !2298, arg: 1, file: !2033, type: !1213, flags: 64)
!2312 = !DILocation(line: 84, column: 1, scope: !2300)
!2313 = !{  }
!2314 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN8T3D_fix3D0Ev", type: !2296, spFlags: 8, unit: !10)
!2315 = !DILocation(scope: !2314)
!2316 = !DILexicalBlock(file: !3, scope: !2314, line: 1, column: 1)
!2317 = !DILocation(scope: !2316)
!2318 = !DILexicalBlock(file: !3, scope: !2316, line: 1, column: 1)
!2319 = !DILocation(scope: !2318)
!2320 = !DILexicalBlock(file: !3, scope: !2318, line: 1, column: 1)
!2321 = !DILocation(scope: !2320)
!2322 = !DILexicalBlock(file: !3, scope: !2320, line: 1, column: 1)
!2323 = !DILocation(scope: !2322)
!2324 = !DILexicalBlock(file: !3, scope: !2316, line: 1, column: 1)
!2325 = !DILocation(scope: !2324)
!2326 = !DILexicalBlock(file: !3, scope: !2324, line: 1, column: 1)
!2327 = !DILocation(scope: !2326)
!2328 = !DILexicalBlock(file: !3, scope: !2326, line: 1, column: 1)
!2329 = !DILocation(scope: !2328)
!2330 = !DILocation(line: 84, column: 1, scope: !2316)
!2331 = !{  }
!2332 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN8T3D_fix3D2Ev", type: !2296, spFlags: 8, unit: !10)
!2333 = !DILocation(scope: !2332)
!2334 = !DILexicalBlock(file: !3, scope: !2332, line: 1, column: 1)
!2335 = !DILocation(scope: !2334)
!2336 = !DILexicalBlock(file: !3, scope: !2334, line: 1, column: 1)
!2337 = !DILocation(scope: !2336)
!2338 = !DILexicalBlock(file: !3, scope: !2336, line: 1, column: 1)
!2339 = !DILocation(scope: !2338)
!2340 = !DILexicalBlock(file: !3, scope: !2338, line: 1, column: 1)
!2341 = !DILocation(scope: !2340)
!2342 = !DILexicalBlock(file: !3, scope: !2334, line: 1, column: 1)
!2343 = !DILocation(scope: !2342)
!2344 = !DILexicalBlock(file: !3, scope: !2342, line: 1, column: 1)
!2345 = !DILocation(scope: !2344)
!2346 = !DILexicalBlock(file: !3, scope: !2344, line: 1, column: 1)
!2347 = !DILocation(scope: !2346)
!2348 = !DILocation(line: 84, column: 1, scope: !2334)
!2349 = !{ null, !1242, !1757, !731, !731, !731 }
!2350 = !DISubroutineType(types: !2349)
!2351 = !{ !2365, !2367, !2369, !2371, !2373 }
!2352 = distinct !DISubprogram(file: !3, scope: !1236, name: "Tinty_f2D", line: 225, type: !2350, spFlags: 8, unit: !10, scopeLine: 225)
!2353 = !DILocation(scope: !2352)
!2354 = !DILexicalBlock(file: !3, scope: !2352, line: 225, column: 1)
!2355 = !DILocation(scope: !2354)
!2356 = !DILexicalBlock(file: !3, scope: !2354, line: 1, column: 1)
!2357 = !DILocation(scope: !2356)
!2358 = !DILexicalBlock(file: !3, scope: !2356, line: 1, column: 1)
!2359 = !DILocation(scope: !2358)
!2360 = !DILexicalBlock(file: !3, scope: !2354, line: 1, column: 1)
!2361 = !DILocation(scope: !2360)
!2362 = !DILexicalBlock(file: !3, scope: !2360, line: 1, column: 1)
!2363 = !DILocation(scope: !2362)
!2364 = !DILocalVariable(scope: !2354, file: !3, type: !1242, flags: 64)
!2365 = !DILocalVariable(scope: !2352, arg: 1, file: !3, type: !1242, flags: 64)
!2366 = !DILocalVariable(scope: !2354, name: "f1", file: !3, type: !1757)
!2367 = !DILocalVariable(scope: !2352, name: "f1", arg: 2, file: !3, type: !1757)
!2368 = !DILocalVariable(scope: !2354, name: "ymin", file: !3, type: !731)
!2369 = !DILocalVariable(scope: !2352, name: "ymin", arg: 3, file: !3, type: !731)
!2370 = !DILocalVariable(scope: !2354, name: "ymax", file: !3, type: !731)
!2371 = !DILocalVariable(scope: !2352, name: "ymax", arg: 4, file: !3, type: !731)
!2372 = !DILocalVariable(scope: !2354, name: "absacc", file: !3, type: !731)
!2373 = !DILocalVariable(scope: !2352, name: "absacc", arg: 5, file: !3, type: !731)
!2374 = !DILocation(line: 225, column: 1, scope: !2354)
!2375 = !{  }
!2376 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN9Tinty_f2DC2ERK11T2DFunctionddd", type: !2350, spFlags: 8, unit: !10)
!2377 = !DILocation(scope: !2376)
!2378 = !DILexicalBlock(file: !3, scope: !2376, line: 1, column: 1)
!2379 = !DILocation(scope: !2378)
!2380 = !DILexicalBlock(file: !3, scope: !2378, line: 1, column: 1)
!2381 = !DILocation(scope: !2380)
!2382 = !DILexicalBlock(file: !3, scope: !2380, line: 1, column: 1)
!2383 = !DILocation(scope: !2382)
!2384 = !DILexicalBlock(file: !3, scope: !2382, line: 1, column: 1)
!2385 = !DILocation(scope: !2384)
!2386 = !DILexicalBlock(file: !3, scope: !2378, line: 1, column: 1)
!2387 = !DILocation(scope: !2386)
!2388 = !DILexicalBlock(file: !3, scope: !2386, line: 1, column: 1)
!2389 = !DILocation(scope: !2388)
!2390 = !DILexicalBlock(file: !3, scope: !2388, line: 1, column: 1)
!2391 = !DILocation(scope: !2390)
!2392 = !DILocation(line: 225, column: 1, scope: !2378)
!2393 = !{ !2487, !2489 }
!2394 = distinct !DISubprogram(file: !3, scope: !10, name: "call", line: 226, type: !1244, spFlags: 8, unit: !10, scopeLine: 226)
!2395 = !DILocation(scope: !2394)
!2396 = !DILexicalBlock(file: !3, scope: !2394, line: 226, column: 1)
!2397 = !DILocation(scope: !2396)
!2398 = !DILexicalBlock(file: !3, scope: !2396, line: 1, column: 1)
!2399 = !DILocation(scope: !2398)
!2400 = !DILexicalBlock(file: !3, scope: !2398, line: 1, column: 1)
!2401 = !DILocation(scope: !2400)
!2402 = !DILexicalBlock(file: !3, scope: !2400, line: 1, column: 1)
!2403 = !DILocation(scope: !2402)
!2404 = !DILexicalBlock(file: !3, scope: !2396, line: 1, column: 1)
!2405 = !DILocation(scope: !2404)
!2406 = !DILexicalBlock(file: !3, scope: !2404, line: 1, column: 1)
!2407 = !DILocation(scope: !2406)
!2408 = !DILexicalBlock(file: !3, scope: !2406, line: 52, column: 1)
!2409 = !DILocation(scope: !2408)
!2410 = !DILexicalBlock(file: !3, scope: !2404, line: 1, column: 1)
!2411 = !DILocation(scope: !2410)
!2412 = !DILexicalBlock(file: !3, scope: !2410, line: 76, column: 1)
!2413 = !DILocation(scope: !2412)
!2414 = !DILexicalBlock(file: !3, scope: !2410, line: 87, column: 1)
!2415 = !DILocation(scope: !2414)
!2416 = !DILexicalBlock(file: !3, scope: !2414, line: 88, column: 1)
!2417 = !DILocation(scope: !2416)
!2418 = !DILexicalBlock(file: !3, scope: !2404, line: 1, column: 1)
!2419 = !DILocation(scope: !2418)
!2420 = !DILexicalBlock(file: !3, scope: !2418, line: 1, column: 1)
!2421 = !DILocation(scope: !2420)
!2422 = !DILexicalBlock(file: !3, scope: !2420, line: 1, column: 1)
!2423 = !DILocation(scope: !2422)
!2424 = !DILexicalBlock(file: !3, scope: !2422, line: 1, column: 1)
!2425 = !DILocation(scope: !2424)
!2426 = !DILexicalBlock(file: !3, scope: !2424, line: 1, column: 1)
!2427 = !DILocation(scope: !2426)
!2428 = !DILexicalBlock(file: !3, scope: !2420, line: 1, column: 1)
!2429 = !DILocation(scope: !2428)
!2430 = !DILexicalBlock(file: !3, scope: !2396, line: 1, column: 1)
!2431 = !DILocation(scope: !2430)
!2432 = !DILexicalBlock(file: !3, scope: !2430, line: 1, column: 1)
!2433 = !DILocation(scope: !2432)
!2434 = !DILexicalBlock(file: !3, scope: !2432, line: 1, column: 1)
!2435 = !DILocation(scope: !2434)
!2436 = !DILexicalBlock(file: !3, scope: !2396, line: 1, column: 1)
!2437 = !DILocation(scope: !2436)
!2438 = !DILexicalBlock(file: !3, scope: !2436, line: 1, column: 1)
!2439 = !DILocation(scope: !2438)
!2440 = !DILexicalBlock(file: !3, scope: !2438, line: 1, column: 1)
!2441 = !DILocation(scope: !2440)
!2442 = !DILexicalBlock(file: !3, scope: !2396, line: 1, column: 1)
!2443 = !DILocation(scope: !2442)
!2444 = !DILexicalBlock(file: !3, scope: !2442, line: 1, column: 1)
!2445 = !DILocation(scope: !2444)
!2446 = !DILexicalBlock(file: !3, scope: !2444, line: 1, column: 1)
!2447 = !DILocation(scope: !2446)
!2448 = !DILexicalBlock(file: !3, scope: !2396, line: 1, column: 1)
!2449 = !DILocation(scope: !2448)
!2450 = !DILexicalBlock(file: !3, scope: !2448, line: 1, column: 1)
!2451 = !DILocation(scope: !2450)
!2452 = !DILexicalBlock(file: !3, scope: !2450, line: 52, column: 1)
!2453 = !DILocation(scope: !2452)
!2454 = !DILexicalBlock(file: !3, scope: !2448, line: 1, column: 1)
!2455 = !DILocation(scope: !2454)
!2456 = !DILexicalBlock(file: !3, scope: !2454, line: 76, column: 1)
!2457 = !DILocation(scope: !2456)
!2458 = !DILexicalBlock(file: !3, scope: !2454, line: 87, column: 1)
!2459 = !DILocation(scope: !2458)
!2460 = !DILexicalBlock(file: !3, scope: !2458, line: 88, column: 1)
!2461 = !DILocation(scope: !2460)
!2462 = !DILexicalBlock(file: !3, scope: !2448, line: 1, column: 1)
!2463 = !DILocation(scope: !2462)
!2464 = !DILexicalBlock(file: !3, scope: !2462, line: 1, column: 1)
!2465 = !DILocation(scope: !2464)
!2466 = !DILexicalBlock(file: !3, scope: !2464, line: 1, column: 1)
!2467 = !DILocation(scope: !2466)
!2468 = !DILexicalBlock(file: !3, scope: !2466, line: 1, column: 1)
!2469 = !DILocation(scope: !2468)
!2470 = !DILexicalBlock(file: !3, scope: !2468, line: 1, column: 1)
!2471 = !DILocation(scope: !2470)
!2472 = !DILexicalBlock(file: !3, scope: !2464, line: 1, column: 1)
!2473 = !DILocation(scope: !2472)
!2474 = !DILexicalBlock(file: !3, scope: !2396, line: 1, column: 1)
!2475 = !DILocation(scope: !2474)
!2476 = !DILexicalBlock(file: !3, scope: !2474, line: 1, column: 1)
!2477 = !DILocation(scope: !2476)
!2478 = !DILexicalBlock(file: !3, scope: !2476, line: 1, column: 1)
!2479 = !DILocation(scope: !2478)
!2480 = !DILexicalBlock(file: !3, scope: !2396, line: 1, column: 1)
!2481 = !DILocation(scope: !2480)
!2482 = !DILexicalBlock(file: !3, scope: !2480, line: 1, column: 1)
!2483 = !DILocation(scope: !2482)
!2484 = !DILexicalBlock(file: !3, scope: !2482, line: 1, column: 1)
!2485 = !DILocation(scope: !2484)
!2486 = !DILocalVariable(scope: !2396, file: !3, type: !1242, flags: 64)
!2487 = !DILocalVariable(scope: !2394, arg: 1, file: !3, type: !1242, flags: 64)
!2488 = !DILocalVariable(scope: !2396, name: "x", file: !3, type: !731)
!2489 = !DILocalVariable(scope: !2394, name: "x", arg: 2, file: !3, type: !731)
!2490 = !DILocation(line: 226, column: 1, scope: !2396)
!2491 = !DILocation(line: 40, column: 1, scope: !2396)
!2492 = !DILocation(line: 184, column: 1, scope: !2396)
!2493 = !DILocation(line: 188, column: 1, scope: !2396)
!2494 = !DILocation(line: 189, column: 1, scope: !2396)
!2495 = !DILocation(line: 190, column: 1, scope: !2396)
!2496 = !DILocation(line: 191, column: 1, scope: !2396)
!2497 = !DILocation(line: 192, column: 1, scope: !2396)
!2498 = !DILocation(line: 194, column: 1, scope: !2396)
!2499 = !DILocation(line: 197, column: 1, scope: !2396)
!2500 = !DILocation(line: 198, column: 1, scope: !2396)
!2501 = !DILocation(line: 199, column: 1, scope: !2396)
!2502 = !DILocation(line: 200, column: 1, scope: !2396)
!2503 = !DILocation(line: 201, column: 1, scope: !2396)
!2504 = !DILocation(line: 202, column: 1, scope: !2396)
!2505 = !DILocation(line: 204, column: 1, scope: !2396)
!2506 = !DILocation(line: 206, column: 1, scope: !2396)
!2507 = !DILocation(line: 210, column: 1, scope: !2396)
!2508 = !DILocation(line: 211, column: 1, scope: !2396)
!2509 = !DILocation(line: 212, column: 1, scope: !2396)
!2510 = !DILocation(line: 215, column: 1, scope: !2396)
!2511 = !DILocation(line: 42, column: 1, scope: !2396)
!2512 = !{ null, !1242 }
!2513 = !DISubroutineType(types: !2512)
!2514 = !{ !2528 }
!2515 = distinct !DISubprogram(file: !3, scope: !10, name: "~Tinty_f2D", type: !2513, spFlags: 8, unit: !10)
!2516 = !DILocation(scope: !2515)
!2517 = !DILexicalBlock(file: !3, scope: !2515, line: 1, column: 1)
!2518 = !DILocation(scope: !2517)
!2519 = !DILexicalBlock(file: !3, scope: !2517, line: 1, column: 1)
!2520 = !DILocation(scope: !2519)
!2521 = !DILexicalBlock(file: !3, scope: !2519, line: 1, column: 1)
!2522 = !DILocation(scope: !2521)
!2523 = !DILexicalBlock(file: !3, scope: !2517, line: 1, column: 1)
!2524 = !DILocation(scope: !2523)
!2525 = !DILexicalBlock(file: !3, scope: !2523, line: 1, column: 1)
!2526 = !DILocation(scope: !2525)
!2527 = !DILocalVariable(scope: !2517, file: !3, type: !1242, flags: 64)
!2528 = !DILocalVariable(scope: !2515, arg: 1, file: !3, type: !1242, flags: 64)
!2529 = !DILocation(line: 198, column: 1, scope: !2517)
!2530 = !DILocation(line: 29, column: 1, scope: !2517)
!2531 = !{  }
!2532 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN9Tinty_f2DD0Ev", type: !2513, spFlags: 8, unit: !10)
!2533 = !DILocation(scope: !2532)
!2534 = !DILexicalBlock(file: !3, scope: !2532, line: 1, column: 1)
!2535 = !DILocation(scope: !2534)
!2536 = !DILexicalBlock(file: !3, scope: !2534, line: 1, column: 1)
!2537 = !DILocation(scope: !2536)
!2538 = !DILexicalBlock(file: !3, scope: !2536, line: 1, column: 1)
!2539 = !DILocation(scope: !2538)
!2540 = !DILexicalBlock(file: !3, scope: !2538, line: 1, column: 1)
!2541 = !DILocation(scope: !2540)
!2542 = !DILexicalBlock(file: !3, scope: !2534, line: 1, column: 1)
!2543 = !DILocation(scope: !2542)
!2544 = !DILexicalBlock(file: !3, scope: !2542, line: 1, column: 1)
!2545 = !DILocation(scope: !2544)
!2546 = !DILexicalBlock(file: !3, scope: !2544, line: 1, column: 1)
!2547 = !DILocation(scope: !2546)
!2548 = !DILocation(line: 29, column: 1, scope: !2534)
!2549 = !{  }
!2550 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN9Tinty_f2DD2Ev", type: !2513, spFlags: 8, unit: !10)
!2551 = !DILocation(scope: !2550)
!2552 = !DILexicalBlock(file: !3, scope: !2550, line: 1, column: 1)
!2553 = !DILocation(scope: !2552)
!2554 = !DILexicalBlock(file: !3, scope: !2552, line: 1, column: 1)
!2555 = !DILocation(scope: !2554)
!2556 = !DILexicalBlock(file: !3, scope: !2554, line: 1, column: 1)
!2557 = !DILocation(scope: !2556)
!2558 = !DILexicalBlock(file: !3, scope: !2556, line: 1, column: 1)
!2559 = !DILocation(scope: !2558)
!2560 = !DILexicalBlock(file: !3, scope: !2552, line: 1, column: 1)
!2561 = !DILocation(scope: !2560)
!2562 = !DILexicalBlock(file: !3, scope: !2560, line: 1, column: 1)
!2563 = !DILocation(scope: !2562)
!2564 = !DILexicalBlock(file: !3, scope: !2562, line: 1, column: 1)
!2565 = !DILocation(scope: !2564)
!2566 = !DILocation(line: 29, column: 1, scope: !2552)
!2567 = !{ null, !1255, !1897, !731, !731, !731, !731, !731 }
!2568 = !DISubroutineType(types: !2567)
!2569 = !{ !2583, !2585, !2587, !2589, !2591, !2593, !2595 }
!2570 = distinct !DISubprogram(file: !3, scope: !1247, name: "Tintxy_f3D", line: 243, type: !2568, spFlags: 8, unit: !10, scopeLine: 243)
!2571 = !DILocation(scope: !2570)
!2572 = !DILexicalBlock(file: !3, scope: !2570, line: 243, column: 1)
!2573 = !DILocation(scope: !2572)
!2574 = !DILexicalBlock(file: !3, scope: !2572, line: 1, column: 1)
!2575 = !DILocation(scope: !2574)
!2576 = !DILexicalBlock(file: !3, scope: !2574, line: 1, column: 1)
!2577 = !DILocation(scope: !2576)
!2578 = !DILexicalBlock(file: !3, scope: !2572, line: 1, column: 1)
!2579 = !DILocation(scope: !2578)
!2580 = !DILexicalBlock(file: !3, scope: !2578, line: 1, column: 1)
!2581 = !DILocation(scope: !2580)
!2582 = !DILocalVariable(scope: !2572, file: !3, type: !1255, flags: 64)
!2583 = !DILocalVariable(scope: !2570, arg: 1, file: !3, type: !1255, flags: 64)
!2584 = !DILocalVariable(scope: !2572, name: "f1", file: !3, type: !1897)
!2585 = !DILocalVariable(scope: !2570, name: "f1", arg: 2, file: !3, type: !1897)
!2586 = !DILocalVariable(scope: !2572, name: "xmin", file: !3, type: !731)
!2587 = !DILocalVariable(scope: !2570, name: "xmin", arg: 3, file: !3, type: !731)
!2588 = !DILocalVariable(scope: !2572, name: "xmax", file: !3, type: !731)
!2589 = !DILocalVariable(scope: !2570, name: "xmax", arg: 4, file: !3, type: !731)
!2590 = !DILocalVariable(scope: !2572, name: "ymin", file: !3, type: !731)
!2591 = !DILocalVariable(scope: !2570, name: "ymin", arg: 5, file: !3, type: !731)
!2592 = !DILocalVariable(scope: !2572, name: "ymax", file: !3, type: !731)
!2593 = !DILocalVariable(scope: !2570, name: "ymax", arg: 6, file: !3, type: !731)
!2594 = !DILocalVariable(scope: !2572, name: "absacc", file: !3, type: !731)
!2595 = !DILocalVariable(scope: !2570, name: "absacc", arg: 7, file: !3, type: !731)
!2596 = !DILocation(line: 243, column: 1, scope: !2572)
!2597 = !{  }
!2598 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN10Tintxy_f3DC2ERK11T3DFunctionddddd", type: !2568, spFlags: 8, unit: !10)
!2599 = !DILocation(scope: !2598)
!2600 = !DILexicalBlock(file: !3, scope: !2598, line: 1, column: 1)
!2601 = !DILocation(scope: !2600)
!2602 = !DILexicalBlock(file: !3, scope: !2600, line: 1, column: 1)
!2603 = !DILocation(scope: !2602)
!2604 = !DILexicalBlock(file: !3, scope: !2602, line: 1, column: 1)
!2605 = !DILocation(scope: !2604)
!2606 = !DILexicalBlock(file: !3, scope: !2604, line: 1, column: 1)
!2607 = !DILocation(scope: !2606)
!2608 = !DILexicalBlock(file: !3, scope: !2600, line: 1, column: 1)
!2609 = !DILocation(scope: !2608)
!2610 = !DILexicalBlock(file: !3, scope: !2608, line: 1, column: 1)
!2611 = !DILocation(scope: !2610)
!2612 = !DILexicalBlock(file: !3, scope: !2610, line: 1, column: 1)
!2613 = !DILocation(scope: !2612)
!2614 = !DILocation(line: 243, column: 1, scope: !2600)
!2615 = !{ !2657, !2659 }
!2616 = distinct !DISubprogram(file: !3, scope: !10, name: "call", line: 244, type: !1257, spFlags: 8, unit: !10, scopeLine: 244)
!2617 = !DILocation(scope: !2616)
!2618 = !DILexicalBlock(file: !3, scope: !2616, line: 244, column: 1)
!2619 = !DILocation(scope: !2618)
!2620 = !DILexicalBlock(file: !3, scope: !2618, line: 1, column: 1)
!2621 = !DILocation(scope: !2620)
!2622 = !DILexicalBlock(file: !3, scope: !2620, line: 1, column: 1)
!2623 = !DILocation(scope: !2622)
!2624 = !DILexicalBlock(file: !3, scope: !2622, line: 1, column: 1)
!2625 = !DILocation(scope: !2624)
!2626 = !DILexicalBlock(file: !3, scope: !2618, line: 1, column: 1)
!2627 = !DILocation(scope: !2626)
!2628 = !DILexicalBlock(file: !3, scope: !2626, line: 1, column: 1)
!2629 = !DILocation(scope: !2628)
!2630 = !DILexicalBlock(file: !3, scope: !2628, line: 1, column: 1)
!2631 = !DILocation(scope: !2630)
!2632 = !DILexicalBlock(file: !3, scope: !2618, line: 1, column: 1)
!2633 = !DILocation(scope: !2632)
!2634 = !DILexicalBlock(file: !3, scope: !2632, line: 1, column: 1)
!2635 = !DILocation(scope: !2634)
!2636 = !DILexicalBlock(file: !3, scope: !2634, line: 1, column: 1)
!2637 = !DILocation(scope: !2636)
!2638 = !DILexicalBlock(file: !3, scope: !2618, line: 1, column: 1)
!2639 = !DILocation(scope: !2638)
!2640 = !DILexicalBlock(file: !3, scope: !2638, line: 1, column: 1)
!2641 = !DILocation(scope: !2640)
!2642 = !DILexicalBlock(file: !3, scope: !2640, line: 1, column: 1)
!2643 = !DILocation(scope: !2642)
!2644 = !DILexicalBlock(file: !3, scope: !2618, line: 1, column: 1)
!2645 = !DILocation(scope: !2644)
!2646 = !DILexicalBlock(file: !3, scope: !2644, line: 1, column: 1)
!2647 = !DILocation(scope: !2646)
!2648 = !DILexicalBlock(file: !3, scope: !2646, line: 1, column: 1)
!2649 = !DILocation(scope: !2648)
!2650 = !DILexicalBlock(file: !3, scope: !2618, line: 1, column: 1)
!2651 = !DILocation(scope: !2650)
!2652 = !DILexicalBlock(file: !3, scope: !2650, line: 1, column: 1)
!2653 = !DILocation(scope: !2652)
!2654 = !DILexicalBlock(file: !3, scope: !2652, line: 1, column: 1)
!2655 = !DILocation(scope: !2654)
!2656 = !DILocalVariable(scope: !2618, file: !3, type: !1255, flags: 64)
!2657 = !DILocalVariable(scope: !2616, arg: 1, file: !3, type: !1255, flags: 64)
!2658 = !DILocalVariable(scope: !2618, name: "z", file: !3, type: !731)
!2659 = !DILocalVariable(scope: !2616, name: "z", arg: 2, file: !3, type: !731)
!2660 = !DILocation(line: 244, column: 1, scope: !2618)
!2661 = !DILocation(line: 82, column: 1, scope: !2618)
!2662 = !DILocation(line: 84, column: 1, scope: !2618)
!2663 = !{ null, !1255 }
!2664 = !DISubroutineType(types: !2663)
!2665 = !{ !2679 }
!2666 = distinct !DISubprogram(file: !3, scope: !10, name: "~Tintxy_f3D", type: !2664, spFlags: 8, unit: !10)
!2667 = !DILocation(scope: !2666)
!2668 = !DILexicalBlock(file: !3, scope: !2666, line: 1, column: 1)
!2669 = !DILocation(scope: !2668)
!2670 = !DILexicalBlock(file: !3, scope: !2668, line: 1, column: 1)
!2671 = !DILocation(scope: !2670)
!2672 = !DILexicalBlock(file: !3, scope: !2670, line: 1, column: 1)
!2673 = !DILocation(scope: !2672)
!2674 = !DILexicalBlock(file: !3, scope: !2668, line: 1, column: 1)
!2675 = !DILocation(scope: !2674)
!2676 = !DILexicalBlock(file: !3, scope: !2674, line: 1, column: 1)
!2677 = !DILocation(scope: !2676)
!2678 = !DILocalVariable(scope: !2668, file: !3, type: !1255, flags: 64)
!2679 = !DILocalVariable(scope: !2666, arg: 1, file: !3, type: !1255, flags: 64)
!2680 = !DILocation(line: 244, column: 1, scope: !2668)
!2681 = !DILocation(line: 29, column: 1, scope: !2668)
!2682 = !{  }
!2683 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN10Tintxy_f3DD0Ev", type: !2664, spFlags: 8, unit: !10)
!2684 = !DILocation(scope: !2683)
!2685 = !DILexicalBlock(file: !3, scope: !2683, line: 1, column: 1)
!2686 = !DILocation(scope: !2685)
!2687 = !DILexicalBlock(file: !3, scope: !2685, line: 1, column: 1)
!2688 = !DILocation(scope: !2687)
!2689 = !DILexicalBlock(file: !3, scope: !2687, line: 1, column: 1)
!2690 = !DILocation(scope: !2689)
!2691 = !DILexicalBlock(file: !3, scope: !2689, line: 1, column: 1)
!2692 = !DILocation(scope: !2691)
!2693 = !DILexicalBlock(file: !3, scope: !2685, line: 1, column: 1)
!2694 = !DILocation(scope: !2693)
!2695 = !DILexicalBlock(file: !3, scope: !2693, line: 1, column: 1)
!2696 = !DILocation(scope: !2695)
!2697 = !DILexicalBlock(file: !3, scope: !2695, line: 1, column: 1)
!2698 = !DILocation(scope: !2697)
!2699 = !DILocation(line: 29, column: 1, scope: !2685)
!2700 = !{  }
!2701 = distinct !DISubprogram(file: !3, scope: !10, name: "_ZN10Tintxy_f3DD2Ev", type: !2664, spFlags: 8, unit: !10)
!2702 = !DILocation(scope: !2701)
!2703 = !DILexicalBlock(file: !3, scope: !2701, line: 1, column: 1)
!2704 = !DILocation(scope: !2703)
!2705 = !DILexicalBlock(file: !3, scope: !2703, line: 1, column: 1)
!2706 = !DILocation(scope: !2705)
!2707 = !DILexicalBlock(file: !3, scope: !2705, line: 1, column: 1)
!2708 = !DILocation(scope: !2707)
!2709 = !DILexicalBlock(file: !3, scope: !2707, line: 1, column: 1)
!2710 = !DILocation(scope: !2709)
!2711 = !DILexicalBlock(file: !3, scope: !2703, line: 1, column: 1)
!2712 = !DILocation(scope: !2711)
!2713 = !DILexicalBlock(file: !3, scope: !2711, line: 1, column: 1)
!2714 = !DILocation(scope: !2713)
!2715 = !DILexicalBlock(file: !3, scope: !2713, line: 1, column: 1)
!2716 = !DILocation(scope: !2715)
!2717 = !DILocation(line: 29, column: 1, scope: !2703)
!2718 = !DIFile(filename: "/usr/include/c++/9/bits/char_traits.h", directory: "/home/talgat/vlasiator")
; !2719 = !DIFile(tag: DW_TAG_file_type, pair: !2718)
!2719 = !{ i32 41, !2718 }
!2720 = !{ !30, !44 }
!2721 = !DISubroutineType(types: !2720)
!2722 = !{ !2728 }
!2723 = distinct !DISubprogram(file: !2718, scope: !445, name: "length", line: 330, type: !2721, spFlags: 8, unit: !10, scopeLine: 330)
!2724 = !DILocation(scope: !2723)
!2725 = !DILexicalBlock(file: !2718, scope: !2723, line: 330, column: 1)
!2726 = !DILocation(scope: !2725)
!2727 = !DILocalVariable(scope: !2725, name: "__s", file: !2718, type: !44)
!2728 = !DILocalVariable(scope: !2723, name: "__s", arg: 1, file: !2718, type: !44)
!2729 = !DILocation(line: 335, column: 1, scope: !2725)
!2730 = !DILocation(line: 336, column: 1, scope: !2725)
!2731 = !DIFile(filename: "/usr/include/c++/9/bits/basic_ios.h", directory: "/home/talgat/vlasiator")
; !2732 = !DIFile(tag: DW_TAG_file_type, pair: !2731)
!2732 = !{ i32 41, !2731 }
!2733 = !DIDerivedType(tag: DW_TAG_pointer_type, size: 64, align: 64, baseType: !118)
!2734 = !{ !162, !2733 }
!2735 = !DISubroutineType(types: !2734)
!2736 = !{ !2742 }
!2737 = distinct !DISubprogram(file: !2731, scope: !118, name: "rdstate", line: 138, type: !2735, spFlags: 8, unit: !10, scopeLine: 138)
!2738 = !DILocation(scope: !2737)
!2739 = !DILexicalBlock(file: !2731, scope: !2737, line: 138, column: 1)
!2740 = !DILocation(scope: !2739)
!2741 = !DILocalVariable(scope: !2739, file: !2731, type: !2733, flags: 64)
!2742 = !DILocalVariable(scope: !2737, arg: 1, file: !2731, type: !2733, flags: 64)
!2743 = !DILocation(line: 138, column: 1, scope: !2739)
!2744 = !{ null, !2733, !162 }
!2745 = !DISubroutineType(types: !2744)
!2746 = !{ !2760, !2762 }
!2747 = distinct !DISubprogram(file: !2731, scope: !118, name: "setstate", line: 158, type: !2745, spFlags: 8, unit: !10, scopeLine: 158)
!2748 = !DILocation(scope: !2747)
!2749 = !DILexicalBlock(file: !2731, scope: !2747, line: 158, column: 1)
!2750 = !DILocation(scope: !2749)
!2751 = !DILexicalBlock(file: !2731, scope: !2749, line: 1, column: 1)
!2752 = !DILocation(scope: !2751)
!2753 = !DILexicalBlock(file: !2731, scope: !2751, line: 1, column: 1)
!2754 = !DILocation(scope: !2753)
!2755 = !DILexicalBlock(file: !2731, scope: !2749, line: 1, column: 1)
!2756 = !DILocation(scope: !2755)
!2757 = !DILexicalBlock(file: !2731, scope: !2755, line: 1, column: 1)
!2758 = !DILocation(scope: !2757)
!2759 = !DILocalVariable(scope: !2749, file: !2731, type: !2733, flags: 64)
!2760 = !DILocalVariable(scope: !2747, arg: 1, file: !2731, type: !2733, flags: 64)
!2761 = !DILocalVariable(scope: !2749, name: "__state", file: !2731, type: !162)
!2762 = !DILocalVariable(scope: !2747, name: "__state", arg: 2, file: !2731, type: !162)
!2763 = !DILocation(line: 158, column: 1, scope: !2749)
!2764 = !DILocation(line: 170, column: 1, scope: !2749)
!2765 = !DIFile(filename: "/usr/include/c++/9/bits/ios_base.h", directory: "/home/talgat/vlasiator")
; !2766 = !DIFile(tag: DW_TAG_file_type, pair: !2765)
!2766 = !{ i32 41, !2765 }
!2767 = !{ !162, !162, !162 }
!2768 = !DISubroutineType(types: !2767)
!2769 = !{ !2775, !2777 }
!2770 = distinct !DISubprogram(file: !2765, scope: !11, name: "operator|", line: 170, type: !2768, spFlags: 8, unit: !10, scopeLine: 170)
!2771 = !DILocation(scope: !2770)
!2772 = !DILexicalBlock(file: !2765, scope: !2770, line: 170, column: 1)
!2773 = !DILocation(scope: !2772)
!2774 = !DILocalVariable(scope: !2772, name: "__a", file: !2765, type: !162)
!2775 = !DILocalVariable(scope: !2770, name: "__a", arg: 1, file: !2765, type: !162)
!2776 = !DILocalVariable(scope: !2772, name: "__b", file: !2765, type: !162)
!2777 = !DILocalVariable(scope: !2770, name: "__b", arg: 2, file: !2765, type: !162)
!2778 = !DILocation(line: 170, column: 1, scope: !2772)
!2779 = !DIFile(filename: "/usr/include/c++/9/ostream", directory: "/home/talgat/vlasiator")
; !2780 = !DIFile(tag: DW_TAG_file_type, pair: !2779)
!2780 = !{ i32 41, !2779 }
!2781 = !{ !2783, !2784 }
!2782 = !DICompositeType(tag: DW_TAG_structure_type, file: !2779, name: "_ZSo", size: 2176, align: 64, elements: !2781)
!2783 = !DIDerivedType(tag: DW_TAG_member, file: !2779, scope: !2782, name: "__vptr", size: 64, align: 64, baseType: !126)
!2784 = !DIDerivedType(tag: DW_TAG_member, file: !2779, scope: !2782, name: "__v_St9basic_iosIcSt11char_traitsIcEE", size: 2112, align: 64, offset: 64, baseType: !118)
!2785 = !DIDerivedType(tag: DW_TAG_reference_type, size: 64, align: 64, baseType: !2782)
!2786 = !{ !216, !2785, !44 }
!2787 = !DISubroutineType(types: !2786)
!2788 = !{ !2810, !2812 }
!2789 = distinct !DISubprogram(file: !2779, scope: !11, name: "operator<<", line: 566, type: !2787, spFlags: 8, unit: !10, scopeLine: 566)
!2790 = !DILocation(scope: !2789)
!2791 = !DILexicalBlock(file: !2779, scope: !2789, line: 566, column: 1)
!2792 = !DILocation(scope: !2791)
!2793 = !DILexicalBlock(file: !2779, scope: !2791, line: 1, column: 1)
!2794 = !DILocation(scope: !2793)
!2795 = !DILexicalBlock(file: !2779, scope: !2793, line: 1, column: 1)
!2796 = !DILocation(scope: !2795)
!2797 = !DILexicalBlock(file: !2779, scope: !2795, line: 1, column: 1)
!2798 = !DILocation(scope: !2797)
!2799 = !DILexicalBlock(file: !2779, scope: !2791, line: 1, column: 1)
!2800 = !DILocation(scope: !2799)
!2801 = !DILexicalBlock(file: !2779, scope: !2791, line: 1, column: 1)
!2802 = !DILocation(scope: !2801)
!2803 = !DILexicalBlock(file: !2779, scope: !2801, line: 1, column: 1)
!2804 = !DILocation(scope: !2803)
!2805 = !DILexicalBlock(file: !2779, scope: !2803, line: 1, column: 1)
!2806 = !DILocation(scope: !2805)
!2807 = !DILexicalBlock(file: !2779, scope: !2791, line: 1, column: 1)
!2808 = !DILocation(scope: !2807)
!2809 = !DILocalVariable(scope: !2791, name: "__out", file: !2779, type: !2785)
!2810 = !DILocalVariable(scope: !2789, name: "__out", arg: 1, file: !2779, type: !2785)
!2811 = !DILocalVariable(scope: !2791, name: "__s", file: !2779, type: !44)
!2812 = !DILocalVariable(scope: !2789, name: "__s", arg: 2, file: !2779, type: !44)
!2813 = !DILocation(line: 567, column: 1, scope: !2791)
!2814 = !DILocation(line: 568, column: 1, scope: !2791)
!2815 = !DILocation(line: 158, column: 1, scope: !2791)
!2816 = !DILocation(line: 570, column: 1, scope: !2791)
!2817 = !DILocation(line: 335, column: 1, scope: !2791)
!2818 = !DILocation(line: 336, column: 1, scope: !2791)
!2819 = !DILocation(line: 572, column: 1, scope: !2791)
!2820 = !DILocation(line: 573, column: 1, scope: !2791)
!2821 = !{  }
!2822 = distinct !DISubprogram(file: !3, scope: !10, name: "__sti___25_backgroundfield_quadr_cpp_6e0a9365", type: !352, spFlags: 8, unit: !10)
!2823 = !DILocation(scope: !2822)
!2824 = !DILexicalBlock(file: !3, scope: !2822, line: 1, column: 1)
!2825 = !DILocation(scope: !2824)
!2826 = !DILocation(line: 74, column: 1, scope: !2824)
!2827 = distinct !DIGlobalVariable(scope: !10, name: "__I___25_backgroundfield_quadr_cpp_6e0a9365", file: !3, line: 7685, type: !61, isDefinition: true)
!2828 = !DIGlobalVariableExpression(var: !2827, expr: !1366)
!2829 = !{ !2831 }
!2830 = !DICompositeType(tag: DW_TAG_structure_type, file: !3, name: "_ZNSt8ios_base4InitE", size: 8, align: 8, elements: !2829)
!2831 = !DIDerivedType(tag: DW_TAG_member, file: !3, scope: !2830, size: 8, align: 8, baseType: !378)
!2832 = distinct !DIGlobalVariable(scope: !10, name: "_ZN47_INTERNAL_25_backgroundfield_quadr_cpp_6e0a9365St8__ioinitE", file: !3, type: !2830, isLocal: true, isDefinition: true)
!2833 = !DIGlobalVariableExpression(var: !2832, expr: !1366)
!2834 = distinct !DIGlobalVariable(scope: !10, name: "__dso_handle", file: !3, type: !17)
!2835 = !DIGlobalVariableExpression(var: !2834, expr: !1366)
!2836 = distinct !DIGlobalVariable(scope: !10, name: "_ZTI11T1DFunction", file: !3, type: !629, isDefinition: true)
!2837 = !DIGlobalVariableExpression(var: !2836, expr: !1366)
!2838 = distinct !DIGlobalVariable(scope: !10, name: "_ZTI11T2DFunction", file: !3, type: !629, isDefinition: true)
!2839 = !DIGlobalVariableExpression(var: !2838, expr: !1366)
!2840 = distinct !DIGlobalVariable(scope: !10, name: "_ZTI8T2D_fix1", file: !3, type: !641, isDefinition: true)
!2841 = !DIGlobalVariableExpression(var: !2840, expr: !1366)
!2842 = distinct !DIGlobalVariable(scope: !10, name: "_ZTI8T3D_fix3", file: !3, type: !641, isDefinition: true)
!2843 = !DIGlobalVariableExpression(var: !2842, expr: !1366)
!2844 = distinct !DIGlobalVariable(scope: !10, name: "_ZTI9Tinty_f2D", file: !3, type: !641, isDefinition: true)
!2845 = !DIGlobalVariableExpression(var: !2844, expr: !1366)
!2846 = distinct !DIGlobalVariable(scope: !10, name: "_ZTI10Tintxy_f3D", file: !3, type: !641, isDefinition: true)
!2847 = !DIGlobalVariableExpression(var: !2846, expr: !1366)
!2848 = !DICompositeType(tag: DW_TAG_array_type, align: 64, baseType: !126, elements: !377)
!2849 = distinct !DIGlobalVariable(scope: !10, name: "_ZTVN10__cxxabiv117__class_type_infoE", file: !3, type: !2848)
!2850 = !DIGlobalVariableExpression(var: !2849, expr: !1366)
!2851 = !DISubrange(count: 14)
!2852 = !{ !2851 }
!2853 = !DICompositeType(tag: DW_TAG_array_type, size: 112, align: 8, baseType: !43, elements: !2852)
!2854 = distinct !DIGlobalVariable(scope: !10, name: "_ZTS11T1DFunction", file: !3, type: !2853, isDefinition: true)
!2855 = !DIGlobalVariableExpression(var: !2854, expr: !1366)
!2856 = distinct !DIGlobalVariable(scope: !10, name: "_ZTS11T2DFunction", file: !3, type: !2853, isDefinition: true)
!2857 = !DIGlobalVariableExpression(var: !2856, expr: !1366)
!2858 = distinct !DIGlobalVariable(scope: !10, name: "_ZTVN10__cxxabiv120__si_class_type_infoE", file: !3, type: !2848)
!2859 = !DIGlobalVariableExpression(var: !2858, expr: !1366)
!2860 = !DISubrange(count: 10)
!2861 = !{ !2860 }
!2862 = !DICompositeType(tag: DW_TAG_array_type, size: 80, align: 8, baseType: !43, elements: !2861)
!2863 = distinct !DIGlobalVariable(scope: !10, name: "_ZTS8T2D_fix1", file: !3, type: !2862, isDefinition: true)
!2864 = !DIGlobalVariableExpression(var: !2863, expr: !1366)
!2865 = distinct !DIGlobalVariable(scope: !10, name: "_ZTS8T3D_fix3", file: !3, type: !2862, isDefinition: true)
!2866 = !DIGlobalVariableExpression(var: !2865, expr: !1366)
!2867 = !DISubrange(count: 11)
!2868 = !{ !2867 }
!2869 = !DICompositeType(tag: DW_TAG_array_type, size: 88, align: 8, baseType: !43, elements: !2868)
!2870 = distinct !DIGlobalVariable(scope: !10, name: "_ZTS9Tinty_f2D", file: !3, type: !2869, isDefinition: true)
!2871 = !DIGlobalVariableExpression(var: !2870, expr: !1366)
!2872 = !DICompositeType(tag: DW_TAG_array_type, size: 104, align: 8, baseType: !43, elements: !245)
!2873 = distinct !DIGlobalVariable(scope: !10, name: "_ZTS10Tintxy_f3D", file: !3, type: !2872, isDefinition: true)
!2874 = !DIGlobalVariableExpression(var: !2873, expr: !1366)
