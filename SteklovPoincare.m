classdef SteklovPoincare < handle
    properties
        A
        b
        x
        
        id
        ids

        Ass
        Aii
        Ais
        Asi
        H

        bs
        bi

        xs
        xi

        rhs

        fcn_fact

    end
    methods
        function obj=SteklovPoincare(A,b,x)
            obj.A=A;
            obj.b=b;
            obj.x=x;
        end
        function prepare(obj,ids,id)
            obj.Ass=obj.A(ids,ids);
            obj.bs=obj.b(ids);
            obj.xs=obj.x(ids);
            for i=1:length(id)
                obj.Aii{i}=obj.A(id{i},id{i});
                obj.Ais{i}=obj.A(id{i},ids);
                obj.Asi{i}=obj.A(ids,id{i});
                obj.H{i}=chol(obj.Aii{i});
                obj.bi{i}=obj.b(id{i});
                obj.xi{i}=obj.x(id{i});
            end
        end
        function prepareRhs(obj,id)
            obj.rhs=obj.bs;
            for i=1:length(id)
                temp=obj.Aii{i}\obj.bi{i};
                temp=obj.Asi{i}*temp;
                obj.rhs=obj.rhs-temp;
            end
        end
        function prepareApp(obj)
            obj.fcn_fact=@(x)( app_fact(x,obj.Ass,obj.H,obj.Asi,obj.Ais) );
        end

        function [x_dd,flag,res,iter]=solve(obj,toll,maxit)
            [x_dd,flag,res,iter]=pcg(obj.fcn_fact,obj.rhs,toll,maxit,[],[],xs);
            obj.xs=x_dd;
        end
    end
end
