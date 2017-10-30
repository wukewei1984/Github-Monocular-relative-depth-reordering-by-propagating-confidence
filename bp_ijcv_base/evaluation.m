%ijcv_eva=[];
%our_eva=[];

for ni=1:97
load(['D:\360云盘\学习\当前工作\dataset\test\bp_ijcv_base\our_bp_result\',num2str(ni),'_ourorder.mat'] );
load(['D:\360云盘\学习\当前工作\dataset\test\bp_ijcv_base\ijcv_result\',num2str(ni),'_ijcvorder.mat']);
load(['D:\360云盘\学习\当前工作\dataset\test\bp_ijcv_base\base_97\',num2str(ni),'_base_order.mat']);

state=regionprops(C,'PixelIdxList' );
m=0;n=0;
for i=1:max(max(C))
    if unique(unique((B(state(i).PixelIdxList))))==i
        m=m+1;
    else
        break;
    end
        
  % if unique(unique((B(state(i).PixelIdList)))==i
   %     n=n+1;
  %  end
end

%ijcv_eva=ijcv_eva+double(m)/double(max(max(C)));
%our_eva=our_eva+double(m)/double(max(max(C)));
our(1,ni)=[double(m)/double(max(max(C)))];

end
%ijcv=ijcv_eva/97;
%our=our_eva/97;

