package sim.simulation;

class SPIRSimulatorEventDrivenTest extends SPIRSimulatorTest {
    @Override protected SPIRSimulator.Engine engine() {
        return SPIRSimulator.Engine.EVENT_DRIVEN;
    }
}
