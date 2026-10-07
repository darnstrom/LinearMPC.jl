#define MPC_BNB_SOURCE_NONE -1 // No solution has been returned since the reset
#define MPC_BNB_SOURCE_SEARCH 0 // Branch and bound without relaxation
#define MPC_BNB_SOURCE_CANDIDATE 1 // Candidate of the warm start
#define MPC_BNB_SOURCE_DEFERRED 2 // Resolved deferred binary controls
extern int bnb_source; // Source of the solution of the latest call (MPC_BNB_SOURCE_*)
void mpc_reset_bnb_deferred(void);
int mpc_bnb_deferred(void);
