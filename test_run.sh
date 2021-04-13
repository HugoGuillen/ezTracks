activate eztracks
echo "Region mode"
python eztracks.py check test_config/config.region.ini
python eztracks.py prepare test_config/config.region.ini
python eztracks.py draw test_config/config.region.ini
echo "Transcript mode"
python eztracks.py prepare test_config/config.trans.ini
python eztracks.py draw test_config/config.trans.ini
echo "Transcript-relative mode"
python eztracks.py prepare test_config/config.trans_relative.ini
python eztracks.py draw test_config/config.trans_relative.ini