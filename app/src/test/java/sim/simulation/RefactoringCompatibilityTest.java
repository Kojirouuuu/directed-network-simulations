package sim.simulation;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

class RefactoringCompatibilityTest {
    @Test
    void preservesOriginalGraphTraversalAndSeededTrajectories() throws Exception {
        assertEquals("732b03a1add91624a3315ccae9a5e81eb810a56d7af3f84e47d7e24d4667f9b8", RefactoringFixture.fingerprint());
    }
}
