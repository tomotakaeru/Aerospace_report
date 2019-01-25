// リスト6.1: ロータリーエンコーダでカウントアップ/ダウンする回路
//a11_rotary.v -- Read from Rotary encoder and count up/down
module a11_rotary (CLK, RSTn, ROTA, ROTB, LEDout, SEG7OUT, SEG7COM);
   input CLK;
   input RSTn;
   input ROTA;
   input ROTB;
   output [9:0] LEDout;
   output [6:0] SEG7OUT;
   output [3:0] SEG7COM;
   reg [9:0]     counter;
   reg [1:0]     regrotA;   // rotary encoder phase A rising edge

   always @ (posedge CLK or negedge RSTn)   // 10bit up /down counter
     begin
    if(RSTn == 1'b0)
      begin
         counter <= 0;
         regrotA[1:0] <= 2'b0;
      end
    else
      begin
         regrotA[1] <= regrotA[0];
         regrotA[0] <= ROTA;
         
         if (regrotA[1:0] == 2'b01) // rising edge
           if (ROTB == 0)           // A rise && B == 0 : ClockWise 
         begin
            if(counter == 'd255)
              counter <= 0;
            else 
              counter <= counter + 10'b1;
         end
           else                    // A rise && B == 1 : Counter ClockWise
         begin
            if(counter == 'd0)
              counter <= 255;
            else 
              counter <= counter - 10'b1;
         end                       // else: !if(ROTB == 0)
      end                          // else: !if(RSTn == 1'b0)
     end                           // always @ (posedge CLK or negedge RSTn)
   assign LEDout[9:0] = counter;
   BIN8to7SEG3 binto7seg3 (CLK, RSTn, counter, SEG7OUT, SEG7COM);
endmodule // a11_rotary
