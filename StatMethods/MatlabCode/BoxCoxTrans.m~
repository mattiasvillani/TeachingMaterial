function yTrans = BoxCoxTrans(y,lambda)

if lambda == 0
    yTrans = log(y);
else
    yTrans = (y.^lambda - 1)/lambda;
end