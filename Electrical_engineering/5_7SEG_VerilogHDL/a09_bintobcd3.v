//リスト5.5:  a09_bintobcd3.v - アップカウンタ出力を BIN8toBCD3 モジュールに入力してダイナミック表示
// 8bit bin counter - display to 3x7seg display

module a09_bintobcd3(CLK,RSTn,LEDout,SEG7OUT,SEG7COM);
input CLK;
input RSTn;
output [9:0] LEDout;
output [6:0] SEG7OUT;
output [3:0] SEG7COM;

reg [7:0] counter;
reg [21:0] prescaler;
wire carryout;

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
        else if(counter == 8'b1111_1111)
            counter <= 0;
        else if( carryout == 1'b1)
            counter <= counter + 8'b1;
        else
            counter <= counter;
    end
    assign LEDout[7:0] = counter;
    assign LEDout [9:8] = 2'b0;

BIN8to7SEG3 binto7seg3 (CLK, RSTn, counter, SEG7OUT, SEG7COM);

endmodule
