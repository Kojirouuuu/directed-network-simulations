package sim.simulation;

class SPIRObserverEventDrivenTest extends SPIRObserverTest {
    @Override protected SPIRSimulator.Engine engine() {
        return SPIRSimulator.Engine.EVENT_DRIVEN;
    }
}
