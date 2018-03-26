classdef grower < handle
    properties
        capacity
        length
        baseArray
    end
    methods
        function obj = grower()
            obj.capacity = 1;
            obj.length = 0;
        end
        function add(obj, val)
            persistent baseArray capacity length
            if(obj.length >= obj.capacity)
                newCapacity = obj.capacity * 2;                
                obj.baseArray(obj.capacity+1:newCapacity) = obj.baseArray(1);
                obj.capacity = newCapacity;
            end            
            obj.length = obj.length+1;
            obj.baseArray(obj.length) = val;
        end
        function a = array(obj)
            a = obj.baseArray(1:obj.length);
        end
    end
end