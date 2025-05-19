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
        function [X,Y,Z]=getCooCenters(obj)
            c=cell(3,1);
            for i=1:3
                c{i}=(0:obj.nlv(i)-1)*obj.h(i);
                c{i}=0.5*(c{i}(1:end-1)+c{i}(2:end));
                if length(c{i})==0
                    c{i}=0;
                end
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
        function ne=get_ne(obj)
            nle=obj.nlv-1;
            nle(nle==0)=1;
            ne=prod(nle);
        end
        function conn=get_connectivity(obj)
            id=1:prod(obj.nlv);
            id=reshape(id,obj.nlv);
            if obj.dim==1
                r{1}=id(1:end-1);
                r{1}=r{1}(:)';
                r{2}=id(2:end);
                r{2}=r{2}(:)';
                conn=[r{1};r{2}];
            end
            if obj.dim==2
                r{1}=id(1:end-1,1:end-1);
                r{1}=r{1}(:)';
                r{2}=id(2:end,1:end-1);
                r{2}=r{2}(:)';
                r{3}=id(1:end-1,2:end);
                r{3}=r{3}(:)';
                r{4}=id(2:end,2:end);
                r{4}=r{4}(:)';
                conn=[r{1};r{2};r{3};r{4}];
            end
            if obj.dim==3
                r{1}=id(1:end-1,1:end-1,1:end-1);
                r{2}=id(2:end,1:end-1,1:end-1);
                r{3}=id(1:end-1,2:end,1:end-1);
                r{4}=id(2:end,2:end,1:end-1);
                r{5}=id(1:end-1,1:end-1,2:end);
                r{6}=id(2:end,1:end-1,2:end);
                r{7}=id(1:end-1,2:end,2:end);
                r{8}=id(2:end,2:end,2:end);
                for i=1:length(r)
                    r{i}=r{i}(:)';
                end
                conn=[r{1};r{2};r{3};r{4};r{5};r{6};r{7};r{8}];
            end
        end
    end
end
         
            
            

