
# calcul de gradient
grad <- function(func, x)
{
	grad = rep(0, length(x))

	for (i in 1:length(x))
	{
		xh    = x;
		xh[i] = xh[i] + 1e-10;
		grad[i] = (func(xh) - func(x))/1e-10;
	}
	return (grad)
}

signe <- function(x)
{
	if(x >= 0){
		return( 1);
	}else{
		return(-1);
	}
}

belongs <- function(k, E)
{
	if(length(E) > 0)
	{
		for(i in 1:length(E))
		{
			if(k == E[i]){
				return (TRUE);
			}
		}
	}
	return (FALSE);
}

# Methode des directions realisables

optim_function <- function(x0, func, A, B)
{
	xk = x0;
	k  = 0 ;

	repeat{

		# Determination des contraintes actives

			Ik = c(); # L'indice des contraintes actives
			current = 1;

			for(i in 1:length(A[,1]))
			{
				if( sum(A[i,]*xk) == B[i] )
				{
					Ik[current] = i;
					current = current + 1;
				}
			}

		# Determination de la direction de descente dk


			Ui <- matrix( 0, nrow = 1, ncol = length(A[1, ]));
			Ui <- Ui[-1,];

			for(i in 1:length(A[1, ]))
			{
				v = rep( 0 , length(A[1, ]));
				v[i] = 1;

				Ui = rbind(Ui, v);
				Ui = rbind(Ui,-v);
			}

			Ci = rep( 1, length(Ui[ ,1]))

			if(length(Ik) > 0)
			{
				for(i in 1:length(Ik))
				{
					Ui = rbind(Ui, A[Ik[i], ]);
					Ci[length(Ci) + 1] = 0;
				}
			}

			# Determination de direction initiale
				d0 = rep( 0, length(Ui[1, ]))

				if(length(Ik) > 0)
				{
					index = length(Ui[ ,1])-length(Ik);
					for(i in 1:length(Ik))
					{
						if( Ui[index+i, ] %*% d0 >= Ci[index+i]){
							for(k in 1:length(Ui[index+i, ]))
							{
								if(Ui[index+i, k] != 0)
								{
									d0[k] <- Ci[index+i]/Ui[index+i, k] + signe(Ci[index+i]/Ui[index+i, k]) * 1e-1;
								}
							}
						}
					}
				}

				grad_f_xk = grad(func, xk);

				current_f <- function(x)
				{
					return (sum(grad_f_xk * x));
				}

				gradient_current_f <- function(x)
				{
					return (grad(current_f,x));
				}

				optimum <- constrOptim(d0, current_f, gradient_current_f, ui = -Ui, ci = -Ci, mu = 1e-06);

				dk <- optimum$par;

		# Determination du pas optimal

			# Determination du pas maximale que peut prendre alpha

				alphak_max = 40;

				for(i in 1:length(A[ ,1]))
				{
					if(sum(A[i, ]*dk) > 0 & !belongs(i,Ik))
					{
						value = (B[i] - sum(A[i, ]*xk))/sum(A[i, ]*dk);
						if(alphak_max > value)
						{
							alphak_max <- value;
						}
					}
				}

			# Determination du pas optimal alpha

				phi <- function(alpha)
				{
					return (func(xk + alpha*dk));
				}

				alphak <- optim(alphak_max/2, phi, lower = 0, upper = alphak_max, method =  "L-BFGS-B");
		# Mise a jour de la solution
		xk = xk + (alphak$par)*dk;
		k = k + 1;

		print(xk)
		# break;

		# Condition d'arrete
		if(sum( grad(func,xk)*dk ) > 0)
		{
			break;
		}
	}
}

f <- function (x)
{
	return (1/2*(x[1]^2 + x[2]^2))
}

A = rbind(c( -1,  1),
		  c(  1,  1),
		  c(  0, -1));

B = c( 7, 5,-2)

I = optim_function(c(-2, 3), f, A, B)




