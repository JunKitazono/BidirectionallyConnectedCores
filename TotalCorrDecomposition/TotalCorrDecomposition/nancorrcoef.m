function [ corr,p ] = nancorrcoef(data)
%�����l(nan)�������đ��֌W�����v�Z����
%function [ corr,p ] = nancorrcoef(data)
nonnan=find(~isnan(sum(data,2)));
[corr,p]=corrcoef(data(nonnan,:));