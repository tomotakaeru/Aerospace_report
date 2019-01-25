	module freeproject (CLK, RSTn, ROTA, ROTB, SW, LEDout, RC);
		input CLK;
		input RSTn;
		input ROTA;
		input ROTB;
		input [2:0] SW;	// 2:ignite, 1:back, 0:forward
		output [9:0] LEDout;	// for output check
		output [3:0] RC;	// output signal partial
		reg [3:0] RC;
		reg [9:0] counter;	// rotary pulse
		reg [2:0] CASE;	// (cw/ccw, rotary amount)
		reg [1:0] regrotA;   // rotary phase A rising edge
		reg [11:0] DRIVE;	// output signal overall
		reg [24:0] timer;
		wire carryout;
		integer i,j;
		integer a,b,c;


	always @ (posedge CLK or negedge RSTn)	// counter -> CASE
		begin
		if (RSTn == 1'b0)
			begin
			counter <= 10'b0;
			regrotA[1:0] <= 2'b0; 
			CASE <= 'b000;
			end
		else
			begin
			regrotA[1] <= regrotA[0];
			regrotA[0] <= ROTA;

			if (regrotA[1:0] == 2'b01)
				if (ROTB == 0)	// CW
					begin
					if (counter == 'd255)
						counter <= 0;
					else 
						counter <= counter + 10'b1;
					end
				else	// CCW
					begin
					if (counter == 'd0)
						counter <= 255;
					else
						counter <= counter - 10'b1;
					end
		  
			if ( counter < 'd127 )	// treat as CW control
				begin
				if (counter > 'd10)
					CASE <= 'b111;
				else
					begin
					if (counter > 'd6)
						CASE <= 'b110;
					else
						if (counter > 'd2)
							CASE <= 'b101;
						else
							CASE <= 'b000;
					end
				end
			else	// treat as CCW control
				if (counter > 'd253)
					CASE <= 'b000;
				else
					begin
					if (counter > 'd249)
						CASE <= 'b001;
					else
						if (counter > 'd245)
							CASE <= 'b010;
						else
							CASE <= 'b011;
					end
			end	// else: if (RSTn != 1'b0)
		end

	// assign LEDout[9:0] = CASE;


	always @ (posedge CLK or negedge RSTn)	// (SW[2:0], CASE) -> DRIVE
		begin
		if (RSTn == 1'b0)
			DRIVE <= 12'b0;
		else
			if (SW[2] == 1)	// IGNITE
				DRIVE <= 12'b111111111111;
			else
				begin
					if (SW[0] == SW[1])		// STOP
						DRIVE <= 12'b0;
					else
						begin
						if (SW[0] == 1)	// move FORWARD
							begin
								case (CASE)
									'b000:	DRIVE <= 12'b010101010101;
									'b001:	DRIVE <= 12'b000100010101;
									'b010:	DRIVE <= 12'b000100010001;
									'b011:	DRIVE <= 12'b000100011001;
									'b101:	DRIVE <= 12'b010101000100;
									'b110:	DRIVE <= 12'b010001000100;
									'b111:	DRIVE <= 12'b011001000100;
									default:	DRIVE <= 12'b000000000000;
								endcase
							end
						else	// move BACK
							begin
								case (CASE)
									'b000:	DRIVE <= 12'b101010101010;
									'b001:	DRIVE <= 12'b001000101010;
									'b010:	DRIVE <= 12'b001000100010;
									'b011:	DRIVE <= 12'b001000100110;
									'b101:	DRIVE <= 12'b101010001000;
									'b110:	DRIVE <= 12'b100010001000;
									'b111:	DRIVE <= 12'b100110001000;
									default:	DRIVE <= 12'b000000000000;
								endcase
							end
						end
						end	// else: if(RSTn != 1'b0)
		end

assign LEDout[9:0] = DRIVE[9:0];
//assign LEDout[1] = ROTA;
//assign LEDout[0] = ROTB;


	always @ (posedge CLK or negedge RSTn)	// CLK30MHz, timer==5M -> carryout=1
		begin
			if (RSTn == 1'b0)
				timer <= 25'b0;
			else if (timer == 25'd5000000)
				timer <= 25'b0;
			else
				timer <= timer + 25'b1;
		end

	assign carryout = (timer == 25'd5000000) ? 1'b1 : 1'b0;


	always @ (posedge carryout or negedge RSTn)	// DRIVE -> RemoteControler in every 1/6sec
		begin
			if (RSTn == 1'b0)
				i = 0;
			else
				i = (i+1)%4;	// i=0~3
				j = 4*i;
				a = j+3;
				b = j+2;
				c = j+1;
				RC[3] <= ~DRIVE[a];
				RC[2] <= ~DRIVE[b];
				RC[1] <= ~DRIVE[c];
				RC[0] <= ~DRIVE[j];
		end
		
	endmodule
		

	  
