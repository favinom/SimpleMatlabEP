classdef PointGrid < handle
    properties
        nlv
        nv
        h
        dim
    end
    methods
        function obj=PointGrid(nlv,h)
            if ~isvector(nlv) && ~isempty(nlv)
                error('nv is not a vector');
            end
            if ~isvector(h) && ~isempty(h)
                error('h is not a vector');
            end
            if length(nlv)~=length(h)
                error('nv and h: not compatible size');
            end
            obj.dim=length(nlv);
            if obj.dim==0
                nlv=[1 1 1];
                h=[1 1 1];
            end
            if obj.dim==1    
                nlv(2)=1;
                h(2)=h(1);
                nlv(3)=1;
                h(3)=h(1);
            end
            if obj.dim==2
                nlv(3)=1;
                h(3)=0.5*(h(1)+h(2));
            end
            obj.nlv=nlv;
            obj.h=h;
            obj.nv=prod(obj.nlv);
        end
        function [X,Y,Z]=getCoo(obj)
            c=cell(3,1);
            for i=1:3
                c{i}=(0:obj.nlv(i)-1)*obj.h(i);
            end
            [X,Y,Z]=ndgrid(c{1},c{2},c{3});
        end
        function h=get_h(obj,i)
            c=(0:obj.nlv(i)-1)*obj.h(i);
            h=diff(c);
        end
        function nv=get_nv(obj)
            nv=obj.nv;
        end
    end
end

            
            
            

