/*
 * Automatically Generated from Mathematica.
 * Tue 24 Oct 2023 17:04:44 GMT-04:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
#include<math.h>
/**
 * Copied from Wolfram Mathematica C Definitions file mdefs.hpp
 * Changed marcos to inline functions (Eric Cousineau)
 */
inline double Power(double x, double y) { return pow(x, y); }
inline double Sqrt(double x) { return sqrt(x); }

inline double Abs(double x) { return fabs(x); }

inline double Exp(double x) { return exp(x); }
inline double Log(double x) { return log(x); }

inline double Sin(double x) { return sin(x); }
inline double Cos(double x) { return cos(x); }
inline double Tan(double x) { return tan(x); }

inline double ArcSin(double x) { return asin(x); }
inline double ArcCos(double x) { return acos(x); }
inline double ArcTan(double x) { return atan(x); }

/* update ArcTan function to use atan2 instead. */
inline double ArcTan(double x, double y) { return atan2(y,x); }

inline double Sinh(double x) { return sinh(x); }
inline double Cosh(double x) { return cosh(x); }
inline double Tanh(double x) { return tanh(x); }

const double E	= 2.71828182845904523536029;
const double Pi = 3.14159265358979323846264;
const double Degree = 0.01745329251994329576924;

inline double Sec(double x) { return 1/cos(x); }
inline double Csc(double x) { return 1/sin(x); }

#endif

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t1035;
  double t1022;
  double t1028;
  double t1040;
  double t1045;
  double t1021;
  double t1051;
  double t1060;
  double t1061;
  double t1092;
  double t1093;
  double t1094;
  double t1095;
  double t1096;
  double t1034;
  double t1041;
  double t1043;
  double t1044;
  double t1062;
  double t1070;
  double t1117;
  double t1114;
  double t1115;
  double t1118;
  double t1122;
  double t1123;
  double t1124;
  double t1132;
  double t1133;
  double t1134;
  double t1135;
  double t1136;
  double t1116;
  double t1119;
  double t1120;
  double t1121;
  double t1125;
  double t1126;
  double t1073;
  double t1089;
  double t1090;
  double t1154;
  double t1155;
  double t1156;
  double t1104;
  double t1100;
  double t1101;
  double t1102;
  double t1103;
  double t1105;
  double t1128;
  double t1129;
  double t1130;
  double t1169;
  double t1170;
  double t1171;
  double t1144;
  double t1140;
  double t1141;
  double t1142;
  double t1143;
  double t1145;
  double t1158;
  double t1159;
  double t1160;
  double t1162;
  double t1163;
  double t1165;
  double t1166;
  double t1167;
  double t1173;
  double t1174;
  double t1175;
  double t1177;
  double t1178;
  double t1180;
  double t1181;
  double t1182;
  double t1227;
  double t1228;
  double t1229;
  double t1231;
  double t1232;
  double t1233;
  double t1247;
  double t1248;
  double t1249;
  double t1251;
  double t1252;
  double t1253;
  double t1193;
  double t1194;
  double t1195;
  double t1189;
  double t1190;
  double t1191;
  double t1107;
  double t1108;
  double t1109;
  double t1110;
  double t1097;
  double t1098;
  double t1099;
  double t1199;
  double t1200;
  double t1209;
  double t1210;
  double t1211;
  double t1205;
  double t1206;
  double t1207;
  double t1147;
  double t1148;
  double t1149;
  double t1150;
  double t1137;
  double t1138;
  double t1139;
  double t1215;
  double t1216;
  double t1157;
  double t1172;
  double t1186;
  double t1187;
  double t1188;
  double t1192;
  double t1196;
  double t1197;
  double t1198;
  double t1201;
  double t1202;
  double t1203;
  double t1204;
  double t1208;
  double t1212;
  double t1213;
  double t1214;
  double t1217;
  double t1218;
  double t1221;
  double t1222;
  double t1223;
  double t1224;
  double t1225;
  double t1230;
  double t1234;
  double t1235;
  double t1237;
  double t1238;
  double t1239;
  double t1241;
  double t1242;
  double t1243;
  double t1244;
  double t1245;
  double t1250;
  double t1254;
  double t1255;
  double t1257;
  double t1258;
  double t1259;
  double t1286;
  double t1287;
  double t1288;
  double t1289;
  double t1290;
  double t1291;
  double t1292;
  double t1293;
  double t1220;
  double t1226;
  double t1236;
  double t1240;
  double t1246;
  double t1256;
  double t1260;
  double t1261;
  double t1091;
  double t1106;
  double t1111;
  double t1112;
  double t1266;
  double t1267;
  double t1268;
  double t1269;
  double t1131;
  double t1146;
  double t1151;
  double t1152;
  double t1272;
  double t1273;
  double t1274;
  double t1275;
  t1035 = Cos(var1[3]);
  t1022 = Cos(var1[4]);
  t1028 = Sin(var1[3]);
  t1040 = Sin(var1[4]);
  t1045 = Sin(var1[2]);
  t1021 = Cos(var1[2]);
  t1051 = t1035*t1022;
  t1060 = -1.*t1028*t1040;
  t1061 = t1051 + t1060;
  t1092 = -1.*t1022;
  t1093 = 1. + t1092;
  t1094 = 0.4*t1093;
  t1095 = 0.64*t1022;
  t1096 = t1094 + t1095;
  t1034 = -1.*t1022*t1028;
  t1041 = -1.*t1035*t1040;
  t1043 = t1034 + t1041;
  t1044 = t1021*t1043;
  t1062 = -1.*t1045*t1061;
  t1070 = t1044 + t1062;
  t1117 = Cos(var1[5]);
  t1114 = Cos(var1[6]);
  t1115 = Sin(var1[5]);
  t1118 = Sin(var1[6]);
  t1122 = t1117*t1114;
  t1123 = -1.*t1115*t1118;
  t1124 = t1122 + t1123;
  t1132 = -1.*t1114;
  t1133 = 1. + t1132;
  t1134 = 0.4*t1133;
  t1135 = 0.64*t1114;
  t1136 = t1134 + t1135;
  t1116 = -1.*t1114*t1115;
  t1119 = -1.*t1117*t1118;
  t1120 = t1116 + t1119;
  t1121 = t1021*t1120;
  t1125 = -1.*t1045*t1124;
  t1126 = t1121 + t1125;
  t1073 = -1.*t1035*t1045;
  t1089 = -1.*t1021*t1028;
  t1090 = t1073 + t1089;
  t1154 = t1021*t1035;
  t1155 = -1.*t1045*t1028;
  t1156 = t1154 + t1155;
  t1104 = t1021*t1061;
  t1100 = t1022*t1028;
  t1101 = t1035*t1040;
  t1102 = t1100 + t1101;
  t1103 = -1.*t1045*t1102;
  t1105 = t1103 + t1104;
  t1128 = -1.*t1117*t1045;
  t1129 = -1.*t1021*t1115;
  t1130 = t1128 + t1129;
  t1169 = t1021*t1117;
  t1170 = -1.*t1045*t1115;
  t1171 = t1169 + t1170;
  t1144 = t1021*t1124;
  t1140 = t1114*t1115;
  t1141 = t1117*t1118;
  t1142 = t1140 + t1141;
  t1143 = -1.*t1045*t1142;
  t1145 = t1143 + t1144;
  t1158 = t1035*t1045;
  t1159 = t1021*t1028;
  t1160 = t1158 + t1159;
  t1162 = t1045*t1043;
  t1163 = t1162 + t1104;
  t1165 = t1021*t1102;
  t1166 = t1045*t1061;
  t1167 = t1165 + t1166;
  t1173 = t1117*t1045;
  t1174 = t1021*t1115;
  t1175 = t1173 + t1174;
  t1177 = t1045*t1120;
  t1178 = t1177 + t1144;
  t1180 = t1021*t1142;
  t1181 = t1045*t1124;
  t1182 = t1180 + t1181;
  t1227 = t1096*t1028;
  t1228 = 0.24*t1035*t1040;
  t1229 = t1227 + t1228;
  t1231 = t1035*t1096;
  t1232 = -0.24*t1028*t1040;
  t1233 = t1231 + t1232;
  t1247 = t1136*t1115;
  t1248 = 0.24*t1117*t1118;
  t1249 = t1247 + t1248;
  t1251 = t1117*t1136;
  t1252 = -0.24*t1115*t1118;
  t1253 = t1251 + t1252;
  t1193 = -1.*t1045*t1043;
  t1194 = -1.*t1021*t1061;
  t1195 = t1193 + t1194;
  t1189 = -1.*t1021*t1035;
  t1190 = t1045*t1028;
  t1191 = t1189 + t1190;
  t1107 = t1096*t1022;
  t1108 = Power(t1040,2);
  t1109 = 0.24*t1108;
  t1110 = t1107 + t1109;
  t1097 = t1096*t1040;
  t1098 = -0.24*t1022*t1040;
  t1099 = t1097 + t1098;
  t1199 = -1.*t1021*t1102;
  t1200 = t1199 + t1062;
  t1209 = -1.*t1045*t1120;
  t1210 = -1.*t1021*t1124;
  t1211 = t1209 + t1210;
  t1205 = -1.*t1021*t1117;
  t1206 = t1045*t1115;
  t1207 = t1205 + t1206;
  t1147 = t1136*t1114;
  t1148 = Power(t1118,2);
  t1149 = 0.24*t1148;
  t1150 = t1147 + t1149;
  t1137 = t1136*t1118;
  t1138 = -0.24*t1114*t1118;
  t1139 = t1137 + t1138;
  t1215 = -1.*t1021*t1142;
  t1216 = t1215 + t1125;
  t1157 = 2.*t1090*t1156;
  t1172 = 2.*t1130*t1171;
  t1186 = Power(t1090,2);
  t1187 = t1090*t1160;
  t1188 = Power(t1156,2);
  t1192 = t1156*t1191;
  t1196 = t1195*t1163;
  t1197 = Power(t1105,2);
  t1198 = Power(t1070,2);
  t1201 = t1200*t1167;
  t1202 = Power(t1130,2);
  t1203 = t1130*t1175;
  t1204 = Power(t1171,2);
  t1208 = t1171*t1207;
  t1212 = t1211*t1178;
  t1213 = Power(t1145,2);
  t1214 = Power(t1126,2);
  t1217 = t1216*t1182;
  t1218 = t1186 + t1187 + t1188 + t1192 + t1196 + t1197 + t1198 + t1201 + t1202 + t1203 + t1204 + t1208 + t1212 + t1213 + t1214 + t1217;
  t1221 = Power(t1035,2);
  t1222 = 0.11*t1221;
  t1223 = Power(t1028,2);
  t1224 = 0.11*t1223;
  t1225 = t1222 + t1224;
  t1230 = -1.*t1229*t1061;
  t1234 = -1.*t1043*t1233;
  t1235 = t1230 + t1234;
  t1237 = t1229*t1102;
  t1238 = t1061*t1233;
  t1239 = t1237 + t1238;
  t1241 = Power(t1117,2);
  t1242 = 0.11*t1241;
  t1243 = Power(t1115,2);
  t1244 = 0.11*t1243;
  t1245 = t1242 + t1244;
  t1250 = -1.*t1249*t1124;
  t1254 = -1.*t1120*t1253;
  t1255 = t1250 + t1254;
  t1257 = t1249*t1142;
  t1258 = t1124*t1253;
  t1259 = t1257 + t1258;
  t1286 = -6.72*t1021;
  t1287 = t1191*t1225;
  t1288 = t1200*t1235;
  t1289 = t1195*t1239;
  t1290 = t1207*t1245;
  t1291 = t1216*t1255;
  t1292 = t1211*t1259;
  t1293 = t1286 + t1287 + t1288 + t1289 + t1290 + t1291 + t1292;
  t1220 = -6.72*t1045;
  t1226 = t1090*t1225;
  t1236 = t1105*t1235;
  t1240 = t1070*t1239;
  t1246 = t1130*t1245;
  t1256 = t1145*t1255;
  t1260 = t1126*t1259;
  t1261 = t1220 + t1226 + t1236 + t1240 + t1246 + t1256 + t1260;
  t1091 = 0.11*t1090;
  t1106 = t1099*t1105;
  t1111 = t1110*t1070;
  t1112 = t1091 + t1106 + t1111;
  t1266 = 0.11*t1191;
  t1267 = t1110*t1195;
  t1268 = t1099*t1200;
  t1269 = t1266 + t1267 + t1268;
  t1131 = 0.11*t1130;
  t1146 = t1139*t1145;
  t1151 = t1150*t1126;
  t1152 = t1131 + t1146 + t1151;
  t1272 = 0.11*t1207;
  t1273 = t1150*t1211;
  t1274 = t1139*t1216;
  t1275 = t1272 + t1273 + t1274;
  p_output1[0]=var2[2]*(-0.5*(t1157 + 2.*t1156*t1160 + 2.*t1070*t1163 + 2.*t1105*t1167 + t1172 + 2.*t1171*t1175 + 2.*t1126*t1178 + 2.*t1145*t1182)*var2[0] - 0.5*t1218*var2[1] - 0.5*t1261*var2[2] - 0.5*t1112*var2[3] - 0.12*t1070*var2[4] - 0.5*t1152*var2[5] - 0.12*t1126*var2[6]);
  p_output1[1]=var2[2]*(-0.5*t1218*var2[0] - 0.5*(t1157 + t1172 + 2.*t1090*t1191 + 2.*t1070*t1195 + 2.*t1105*t1200 + 2.*t1130*t1207 + 2.*t1126*t1211 + 2.*t1145*t1216)*var2[1] - 0.5*t1293*var2[2] - 0.5*t1269*var2[3] - 0.12*t1195*var2[4] - 0.5*t1275*var2[5] - 0.12*t1211*var2[6]);
  p_output1[2]=(-0.5*t1261*var2[0] - 0.5*t1293*var2[1])*var2[2];
  p_output1[3]=(-0.5*t1112*var2[0] - 0.5*t1269*var2[1])*var2[2];
  p_output1[4]=(-0.12*t1070*var2[0] - 0.12*t1195*var2[1])*var2[2];
  p_output1[5]=(-0.5*t1152*var2[0] - 0.5*t1275*var2[1])*var2[2];
  p_output1[6]=(-0.12*t1126*var2[0] - 0.12*t1211*var2[1])*var2[2];
}



#ifdef MATLAB_MEX_FILE

#include "mex.h"
/*
 * Main function
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  size_t mrows, ncols;

  double *var1,*var2;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 2)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Two input(s) required (var1,var2).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "Ce2_vec3_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
