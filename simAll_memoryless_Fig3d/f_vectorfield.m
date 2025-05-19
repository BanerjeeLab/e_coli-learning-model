function f = f_vectorfield(t,y,Probl)

if isempty(Probl.param)
    f = feval(Probl.f_fun,t,y) ;
else
    f = feval(Probl.f_fun,t,y,Probl.param) ;
end

end