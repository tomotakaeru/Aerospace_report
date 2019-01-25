// リスト6.4: プッシュスイッチによるカウンタ (完全版)
//a12_pushSW.v -- Push1 - countup Push0 - countdown (Ver.3)

module a12_pushSW (CLK, RSTn, PUSH, LEDout, SEG7OUT, SEG7COM);
   
   input CLK;
   input RSTn;
   input [1:0] PUSH;
   output [9:0] LEDout;
   output [6:0] SEG7OUT;
   output [3:0] SEG7COM;
   reg [9:0]     counter;
   wire     carryout;
   reg [1:0]     regpush1;
   reg [1:0]     regpush0;
   reg [21:0]     prescaler;

   //prescaler for anti-chattering
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
   always @ (posedge CLK or negedge RSTn)
     begin
    if(RSTn == 1'b0)
      begin
         counter <= 0;
         regpush1 [1:0] = 2'b0;
         regpush0 [1:0] = 2'b0;
      end
    else if ( carryout == 1'b1)
      begin
         regpush1[1] <= regpush1[0] ;
         regpush1[0] <= PUSH[1] ;
         regpush0[1] <= regpush0[0] ;
         regpush0[0] <= PUSH[0] ;

         if (regpush1[1:0] == 2'b01) // Push1: countup
           begin
          if(counter == 'd255)
            counter <= 0;
          else 
            counter <= counter + 10'b1;
           end
         if (regpush0[1:0] == 2'b01) // Push0: countdown
           begin
          if(counter == 'd0)
            counter <= 255;
          else 
            counter <= counter - 10'b1;
           end
      end // if ( carryout == 1'b1)
     end                             // always @ (posedge CLK or negedge RSTn)
   
   assign LEDout[9:0] = counter;
   BIN8to7SEG3 binto7seg3 (CLK, RSTn, counter, SEG7OUT, SEG7COM);
endmodule                            // a12_pushSW
