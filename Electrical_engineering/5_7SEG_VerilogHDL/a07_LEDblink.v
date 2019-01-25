//  リスト5.1 10連バー LED を0.1秒間隔で点滅(カウントアップ)する回路

module a07_LEDblink (
    input        CLK,
    input        RSTn,
    output    [9:0]    LEDOUT
);

    reg    [24:0]    prescaler;
    reg    [9:0]    counter;

    wire        carryout;


    always @ (posedge CLK or negedge RSTn) begin
        if(RSTn == 1'b0)
            prescaler <= 25'b0;
        else if(prescaler == 25'd3000000)
            prescaler <= 25'b0;
        else
            prescaler <= prescaler + 25'b1;
    end
    assign    carryout = (prescaler == 25'd3000000) ? 1 : 0;

    always @ (posedge CLK or negedge RSTn) begin
        if(RSTn == 1'b0)
            counter <= 10'b0;
        else if(counter == 10'b11_1111_1111)
            counter <= 10'b0;
        else if( carryout == 1'b1)
            counter <= counter + 10'b1;
        else
            counter <= counter;
    end

    assign    LEDOUT = counter;
endmodule
