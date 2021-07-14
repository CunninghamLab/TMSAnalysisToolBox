function MovingVerticalLineExample

    close all;
    
    hFig = figure;
    xx = -2*pi:0.001:2*pi;
    hAxes = subplot(1,1,1);
    plot(xx,sin(xx),'g');
    
    set(hAxes,'Color','k');
    
    set(hFig,'WindowButtonDownFcn',  @mouseDown);
    set(hFig,'WindowButtonMotionFcn',@mouseMove);
    set(hFig,'WindowButtonUpFcn',    @mouseUp);
    
    randX = max(xx)*rand;
    hVerticalLines = line([randX,randX],get(hAxes,'Ylim'),'Color','red');
    randX = max(xx)*rand;
    hVerticalLines = [hVerticalLines line([randX,randX],get(hAxes,'Ylim'),'Color','blue')];
    
    hLineToDrag = [];
    
    function mouseDown(hObject,~)
        
        % is the mouse down event within the axes?
        if IsCursorInControl(hObject, hAxes)
        
            currentPoint   = get(hAxes,'CurrentPoint');
            xCurrentPoint  = currentPoint(2,1);
            
            for k=1:length(hVerticalLines)
            
                xVertLineCoord = get(hVerticalLines(k),'XData');

                if abs(xCurrentPoint - xVertLineCoord(1)) < 0.1
                    hLineToDrag = hVerticalLines(k);
                    break;
                end
            end
        end
    end

    function mouseUp(~,~)
        
        hLineToDrag = [];
        
    end
    
    function mouseMove(hObject,~)
          
        % is the mouse down event within the axes?
        if ~isempty(hLineToDrag) && IsCursorInControl(hObject, hAxes)
            currentPoint = get(hAxes,'CurrentPoint');
            x            = currentPoint(2,1);
            set(hLineToDrag, 'XData', [x x]);
        end
    end

    function [status] = IsCursorInControl(hCursor, hControl)
        
        status = false;
        
        % get the position of the mouse
        figCurrentPoint = get(hCursor, 'CurrentPoint');
        position      = get(hCursor, 'Position');
        xCursor       = figCurrentPoint(1,1)/position(1,3); % normalize
        yCursor       = figCurrentPoint(1,2)/position(1,4); % normalize

        % get the position of the axes within the GUI
        controlPos = get(hControl,'Position');
        minx    = controlPos(1);
        miny    = controlPos(2);
        maxx    = minx + controlPos(3);
        maxy    = miny + controlPos(4);
            
        % is the mouse down event within the axes?
        if xCursor >= minx && xCursor <= maxx && yCursor >= miny && yCursor <= maxy
            status = true;
        end
    end
end