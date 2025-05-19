###############################################################################################
# Transformation functions - from vinedist [package rvinecopulib] to RVM [package VineCopula] #
###############################################################################################

transformation_to_RVM_all = function(vinecopdist){
  d = vinecopdist$structure$d
  order = vinecopdist$structure$order
  array = vinecopdist$structure$struct_array
  vine = vinecopdist
  p = c((d-1):1)
  
  Matrix = diag(order)
  Family = matrix(0,d,d)
  Par = matrix(0,d,d)
  Par2 = matrix(0,d,d)
  
  for (i in (d-1):1){
    for (j in 1:i){
      ij = order[array[[p[i]]][j]]
      Matrix[(1+i),j]=ij
      
      cop = get_pair_copula(vine,p[i],j)
      Family[(1+i),j] = get_fam_par2(cop)[1]
      Par[(1+i),j] = get_fam_par2(cop)[2]
      Par2[(1+i),j] = get_fam_par2(cop)[3]
    }
  }
  
  RVM = RVineMatrix(Matrix = Matrix, family = Family,
                     par = Par, par2 = Par2)
  
  return(RVM)
}

get_fam_par2 = function(paircop){
  fam = paircop$family
  rot = paircop$rotation
  params = paircop$parameters
  r = 0
  z = 1
  f = 0
  p = 0
  p2 = 0
  
  if (rot == 180){
    r = 1
    z = 1
  }else if (rot == 90){
    r = 3
    z = -1
  }else if (rot == 270){
    r = 2
    z = -1
  }else{
    r = 0
    z = 1
  }
  
  if (fam == "indep"){
    f = 0
    p = 0
    p2 = 0
  }else if (fam == "gaussian"){
    f = 1
    p = params
    p2 = 0
  }else if (fam == "student" || fam == "t"){
    f = 2
    p = params[1]
    p2 = params[2]
  }else if (fam == "clayton"){
    f = 3 + r*10
    p = params*z
    p2 = 0
  }else if (fam == "gumbel"){
    f = 4 + r*10
    p = params*z
    p2 = 0
  }else if (fam == "frank"){
    f = 5
    p = params*z
    p2 = 0
  }else if (fam == "joe"){
    f = 6 + r*10
    p = params*z
    p2 = 0
  }else if (fam == "bb1"){
    f = 7 + r*10
    p = params[1]*z
    p2 = params[2]*z
  }else if (fam == "bb6"){
    f = 8 + r*10
    p = params[1]*z
    p2 = params[2]*z
  }else if (fam == "bb7"){
    f = 9 + r*10
    p = params[1]*z
    p2 = params[2]*z
  }else if (fam == "bb8"){
    f = 10 + r*10
    p = params[1]*z
    p2 = params[2]*z
  }
  
  return(c(f,p,p2))
}
