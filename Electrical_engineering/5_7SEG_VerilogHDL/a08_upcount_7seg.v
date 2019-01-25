// a08_upcount_7seg -  10bit up counter + dynamic 7seg display

module a08_upcount_7seg (CLK, RSTn, LEDout, seg7out, seg7com);

    input CLK;
    input RSTn;
    output [9:0] LEDout;
    output [6:0] seg7out;
    output [3:0] seg7com;

    wire [3:0] seg7com;
    reg [9:0] counter;     // 10bit counter register
    reg [21:0] prescaler;  // prescaler for counter
    wire carryout;          // countup signal for counter
    reg [1:0] counter_7seg;    // 2bit counter for 7seg display
    reg [17:0] prescaler_7seg; // prescaler for 7seg display
    wire carryout_7seg;         // countup signal for 7seg display
    wire [3:0] bin4;  // 4bit binary input to 7seg driver
    
    // 25bit up counter with carry out at 'd3000000
    always @ (posedge CLK or negedge RSTn) begin
        if(RSTn == 1'b0)
            prescaler <= 22'b0;
        else if(prescaler == 22'd3000000)
            prescaler <= 22'b0;
        else
            prescaler <= prescaler + 22'b1;
    end
    assign    carryout = (prescaler == 22'd3000000) ? 1'b1 : 1'b0;

    // 10bit up counter that counts up at carryout = 1;
    always @ (posedge CLK or negedge RSTn) begin
        if(RSTn == 1'b0)
            counter <= 0;
        else if(counter == 10'b11_1111_1111)
            counter <= 0;
        else if( carryout == 1'b1)
            counter <= counter + 10'b1;
        else
            counter <= counter;
    end
    assign LEDout = counter;

// 7segment display part
    // prescaler that countsup x12 times as counter
    always @ (posedge CLK or negedge RSTn) begin
        if(RSTn == 1'b0)
            prescaler_7seg <= 0;
        else if(prescaler_7seg == 'd250000)
            prescaler_7seg <= 18'b0;
        else
            prescaler_7seg <= prescaler_7seg + 18'b1;
    end
    assign carryout_7seg = (prescaler_7seg == 'd250000)? 1'b1 : 1'b0;

    // Display output selector
    always @(posedge CLK or negedge RSTn) begin
        if(RSTn == 1'b0)
            counter_7seg <= 2'b0;
        else if (carryout_7seg == 1'b1)
            if (counter_7seg == 2'b10)
                counter_7seg <= 2'b0;
            else
                counter_7seg <= counter_7seg + 2'b1;
        else
                counter_7seg <= counter_7seg;
    end
    
function [3:0] select7seg;
    input [1:0] counter_7seg;
    input [9:0] counter;
        
    case (counter_7seg)
        default:    select7seg = 4'b1111;
        2'b00:    select7seg = counter [3:0];
        2'b01:    select7seg = counter [7:4];
        2'b10:    begin
                    select7seg[1:0] = counter [9:8];
                    select7seg[3:2] = 2'b0;
                end
    endcase
endfunction

BIN4to7SEG out7seg (select7seg(counter_7seg,counter),seg7out);

function [3:0] select7segcom;
    input [1:0] counter_7seg;
    case (counter_7seg)
        default:    select7segcom = 4'b1111;
        2'b00:    select7segcom = 4'b1110;
        2'b01:    select7segcom = 4'b1101;
        2'b10:    select7segcom = 4'b1011;
    endcase
endfunction
    assign seg7com = select7segcom(counter_7seg);
endmodule
