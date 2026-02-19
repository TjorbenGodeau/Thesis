// filepath: c:/Users/tjorb/Documents/Thesis/rtl/GbSB.sv
// GbSB (Generalized bifurcation Spin-Boson) Hardware Architecture
// Modular design with control interface for external scheduling

module GbSB #(
    parameter int N = 8,           // Number of oscillators
    parameter int DATA_WIDTH = 32, // Fixed-point data width
    parameter int FRAC_WIDTH = 16, // Fractional bits
    parameter bit PER_SPIN = 1'b0  // Individual p_i per spin
) (
    input  logic clk,
    input  logic rst_n,
    
    // Control signals
    input  logic start,
    input  logic [31:0] M,         // Total schedule steps
    output logic done,
    output logic [31:0] step_count,
    
    // Initial conditions
    input  logic init_valid,
    input  logic [N-1:0][DATA_WIDTH-1:0] x0,
    input  logic [N-1:0][DATA_WIDTH-1:0] y0,
    input  logic [N-1:0][DATA_WIDTH-1:0] p0,
    
    // Coupling matrix J (external storage)
    input  logic [N-1:0][N-1:0][DATA_WIDTH-1:0] J,
    
    // Time step control
    input  logic [DATA_WIDTH-1:0] dt,
    input  logic [DATA_WIDTH-1:0] A,  // Nonlinear control constant
    
    // State outputs (synchronous)
    output logic valid_out,
    output logic [N-1:0][DATA_WIDTH-1:0] x,
    output logic [N-1:0][DATA_WIDTH-1:0] y,
    output logic [DATA_WIDTH-1:0] p_global,
    output logic [N-1:0][DATA_WIDTH-1:0] p_individual,
    
    // Results
    output logic [N-1:0] spins,           // ±1 as sign(x)
    output logic [DATA_WIDTH-1:0] energy,
    output logic [DATA_WIDTH-1:0] hamiltonian,
    
    // Status signals
    output logic hit_boundary,            // Inelastic wall collision
    output logic [N-1:0] boundary_mask
);

    // Internal state registers
    logic [N-1:0][DATA_WIDTH-1:0] x_reg, x_next;
    logic [N-1:0][DATA_WIDTH-1:0] y_reg, y_next;
    logic [DATA_WIDTH-1:0] p_global_reg, p_global_next;
    logic [N-1:0][DATA_WIDTH-1:0] p_individual_reg, p_individual_next;
    logic [31:0] m_counter, m_next;
    
    // Step state machine
    typedef enum logic [2:0] {
        IDLE, 
        COMPUTE_Jx,      // Sum J*x
        UPDATE_Y,         // Update y momentum
        UPDATE_X,         // Update position
        CHECK_BOUNDARY,   // Check wall collisions
        UPDATE_P,         // Update bifurcation parameter
        OUTPUT
    } state_t;
    
    state_t state, state_next;

    // =========================================================================
    // Submodule instances
    // =========================================================================
    
    // Matrix-vector multiplication: Jx = J*x
    MatrixVectorMult #(
        .N(N),
        .DATA_WIDTH(DATA_WIDTH),
        .FRAC_WIDTH(FRAC_WIDTH)
    ) jx_mult (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(state == COMPUTE_Jx),
        .J(J),
        .x(x_reg),
        .result(Jx_result),
        .valid_out(jx_valid)
    );
    logic [N-1:0][DATA_WIDTH-1:0] Jx_result;
    logic jx_valid;
    
    // Update Y: y_next = y - dt*(p*x - Jx)
    UpdateMomentum #(
        .N(N),
        .DATA_WIDTH(DATA_WIDTH),
        .FRAC_WIDTH(FRAC_WIDTH)
    ) update_y_block (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(state == UPDATE_Y),
        .y(y_reg),
        .x(x_reg),
        .Jx(Jx_result),
        .p(PER_SPIN ? p_individual_reg : {N{p_global_reg}}),
        .dt(dt),
        .y_next(y_next),
        .valid_out()
    );
    
    // Update X: x_next = x + dt*y_next
    UpdatePosition #(
        .N(N),
        .DATA_WIDTH(DATA_WIDTH),
        .FRAC_WIDTH(FRAC_WIDTH)
    ) update_x_block (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(state == UPDATE_X),
        .x(x_reg),
        .y(y_next),
        .dt(dt),
        .x_next(x_next),
        .valid_out()
    );
    
    // Boundary enforcement: |x| > 1 -> x = sign(x), y = 0
    BoundaryCheck #(
        .N(N),
        .DATA_WIDTH(DATA_WIDTH),
        .FRAC_WIDTH(FRAC_WIDTH)
    ) boundary_check_block (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(state == CHECK_BOUNDARY),
        .x(x_next),
        .y(y_next),
        .x_enforced(x_next),
        .y_enforced(y_next),
        .collision_mask(boundary_mask),
        .any_collision(hit_boundary),
        .valid_out()
    );
    
    // Update P: bifurcation parameter schedule
    UpdateBifurcation #(
        .N(N),
        .DATA_WIDTH(DATA_WIDTH),
        .PER_SPIN(PER_SPIN),
        .FRAC_WIDTH(FRAC_WIDTH)
    ) update_p_block (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(state == UPDATE_P),
        .p_global(p_global_reg),
        .p_individual(p_individual_reg),
        .x(x_reg),
        .m(m_counter),
        .M(M),
        .A(A),
        .p_global_next(p_global_next),
        .p_individual_next(p_individual_next),
        .valid_out()
    );
    
    // Energy computation
    EnergyCompute #(
        .N(N),
        .DATA_WIDTH(DATA_WIDTH),
        .FRAC_WIDTH(FRAC_WIDTH)
    ) energy_block (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(state == OUTPUT),
        .x(x_reg),
        .J(J),
        .energy_out(energy),
        .valid_out()
    );
    
    // Hamiltonian computation
    HamiltonianCompute #(
        .N(N),
        .DATA_WIDTH(DATA_WIDTH),
        .FRAC_WIDTH(FRAC_WIDTH)
    ) hamiltonian_block (
        .clk(clk),
        .rst_n(rst_n),
        .valid_in(state == OUTPUT),
        .x(x_reg),
        .y(y_reg),
        .p_global(p_global_reg),
        .p_individual(p_individual_reg),
        .J(J),
        .hamiltonian_out(hamiltonian),
        .valid_out()
    );

    // =========================================================================
    // State machine and sequential logic
    // =========================================================================
    
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= IDLE;
            x_reg <= '0;
            y_reg <= '0;
            p_global_reg <= '0;
            p_individual_reg <= '0;
            m_counter <= '0;
        end else begin
            state <= state_next;
            x_reg <= x_next;
            y_reg <= y_next;
            p_global_reg <= p_global_next;
            p_individual_reg <= p_individual_next;
            m_counter <= m_next;
        end
    end
    
    always_comb begin
        state_next = state;
        x_next = x_reg;
        y_next = y_reg;
        p_global_next = p_global_reg;
        p_individual_next = p_individual_reg;
        m_next = m_counter;
        done = 1'b0;
        valid_out = 1'b0;
        
        case (state)
            IDLE: begin
                if (start) begin
                    state_next = COMPUTE_Jx;
                    if (init_valid) begin
                        x_next = x0;
                        y_next = y0;
                        p_global_next = p0[0];
                        p_individual_next = p0;
                    end
                    m_next = 32'h0;
                end
            end
            
            COMPUTE_Jx:
                if (jx_valid) state_next = UPDATE_Y;
            
            UPDATE_Y:
                state_next = UPDATE_X;
            
            UPDATE_X:
                state_next = CHECK_BOUNDARY;
            
            CHECK_BOUNDARY:
                state_next = UPDATE_P;
            
            UPDATE_P: begin
                m_next = m_counter + 1;
                if (m_next >= M)
                    state_next = OUTPUT;
                else
                    state_next = COMPUTE_Jx;
            end
            
            OUTPUT: begin
                valid_out = 1'b1;
                done = 1'b1;
                state_next = IDLE;
            end
        endcase
    end
    
    // Output assignments
    assign x = x_reg;
    assign y = y_reg;
    assign p_global = p_global_reg;
    assign p_individual = p_individual_reg;
    assign step_count = m_counter;
    
    // Spin extraction: sign(x)
    generate
        for (genvar i = 0; i < N; i++) begin : gen_spins
            assign spins[i] = (x_reg[i][DATA_WIDTH-1] == 1'b1) ? 1'b0 : 1'b1;
        end
    endgenerate

endmodule