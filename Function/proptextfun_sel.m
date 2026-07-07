function proptextfun_sel(sat,graph,value)

for i=1:sat.node.total_node
    if sat.node.globe(i).ID_item == value
        proptext.node(i).text=text(graph,sat.node.globe(i).node(1),...
            sat.node.globe(i).node(2),sat.node.globe(i).node(3)...
            ,num2str(sat.node.globe(i).ID),'Visible','on','Color','black','BackgroundColor','none','FontWeight','bold');
        hold(graph,'on');
    end
end
end