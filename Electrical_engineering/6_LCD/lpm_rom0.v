// lpm_rom0: LCD module data
module     lpm_rom0 (romaddr, CLK, romq);  // 9bit 32word ROM
	input CLK;
	input [6:0] romaddr;
	output [8:0] romq;

    wire [6:0]  romaddr;  // ROM address
    reg [8:0] romq;     // ROM data

	 always @ (romaddr)
		case (romaddr)
		7'b00000: romq <= 'h038  ;
		7'b00001: romq <= 'h00F  ;
		7'b00010: romq <= 'h001  ;
		7'b00011: romq <= 'h150  ;//from here input ASCII character
		7'b00100: romq <= 'h16C  ;
		7'b00101: romq <= 'h165  ;
		7'b00110: romq <= 'h161  ;
		7'b00111: romq <= 'h173  ;
		7'b01000: romq <= 'h165  ;
		7'b01001: romq <= 'h12C  ;
		7'b01010: romq <= 'h164  ;
		7'b01011: romq <= 'h16F  ;
		7'b01100: romq <= 'h16E  ;
		7'b01101: romq <= 'h127  ;
		7'b01110: romq <= 'h174  ;
		7'b01111: romq <= 'h120  ;
		7'b10000: romq <= 'h173  ;
		7'b10001: romq <= 'h168  ;
		7'b10010: romq <= 'h175  ;
		7'b10011: romq <= 'h0C0  ;//改行
		7'b10100: romq <= 'h174  ;
		7'b10101: romq <= 'h164  ;
		7'b10110: romq <= 'h16F  ;
		7'b10111: romq <= 'h177  ;
		7'b11000: romq <= 'h16E  ;
		7'b11001: romq <= 'h121  ;
		7'b11010: romq <= 'h121  ;
		7'b11011: romq <= 'h121  ;
		7'b11100: romq <= 'h121  ;
		7'b11101: romq <= 'h0FF  ;
		7'b11110: romq <= 'h0  ;
		7'b11111: romq <= 'h0  ;
		endcase
endmodule