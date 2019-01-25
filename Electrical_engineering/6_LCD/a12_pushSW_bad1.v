//リスト6.2: プッシュスイッチによるカウンタ (不完全版1)
//a12_pushSW.v -- Push1 - countup Push0 - countdown (Ver.1)
module a12_pushSW (CLK, RSTn, PUSH, LEDout, SEG7OUT, SEG7COM);
   input CLK;
   input RSTn;
   input [1:0] PUSH;
   output [9:0] LEDout;
   output [6:0] SEG7OUT;
   output [3:0] SEG7COM;
   reg [9:0]     counter;
   always @ (posedge CLK or negedge RSTn)
     begin
    if(RSTn == 1'b0)
      begin
         counter <= 0;
      end
    else
      begin
         if (PUSH[1] == 1) // Push1: countup
           begin
          if(counter == 'd255)
            counter <= 0;
          else 
            counter <= counter + 10'b1;
           end
         if (PUSH[0] == 1) // Push0: countdown
           begin
          if(counter == 'd0)
            counter <= 255;
          else 
            counter <= counter - 10'b1;
           end
      end                  // else: !if(RSTn == 1'b0)
     end                   // always @ (posedge CLK or negedge RSTn)
   assign LEDout[9:0] = counter;
   BIN8to7SEG3 binto7seg3 (CLK, RSTn, counter, SEG7OUT, SEG7COM);
endmodule                  // a12_pushSW
