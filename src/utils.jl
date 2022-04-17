function zoh(A,B,Ts)
  dims = size(B);
  if(length(dims)==1)
	nx = dims[1]
	nu = 1
  else
	nx = dims[1]
	nu = dims[2]
  end
  M = exp([A*Ts  B*Ts;
		   zeros(nu, nx + nu)])
  F = M[1:nx, 1:nx]
  G = M[1:nx, nx+1:nx+nu]
  return F,G
end

